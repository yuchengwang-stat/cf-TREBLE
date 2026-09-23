## cf-TREBLE pipeline -- marker CpG selection and signature estimation ---------
##
## Marker selection.  Nothing here is hardcoded to a particular tree: the node
## counts, the complement columns and the one-past-the-last-complement bound all
## come off the tree bundle, and every threshold comes from config$markers.

beta_IQR <- function(mu, sigma, alpha_cap = 150000) {
  p <- pnorm(mu / sqrt(1 + sigma^2))
  a <- p * (1 - p) / (pbv::pbvnorm(x = mu / sqrt(1 + sigma^2),
                                   y = mu / sqrt(1 + sigma^2),
                                   rho = sigma^2 / (1 + sigma^2)) - p^2) - 1
  a <- pmin(a, alpha_cap)
  qbeta(0.75, p * a, (1 - p) * a) - qbeta(0.25, p * a, (1 - p) * a)
}

## joint node probabilities along the tree
node_joint_p <- function(tree, p) {
  p2 <- p
  walk <- function(j, acc) {
    p2[, j] <<- acc * p2[, j]
    if (length(tree$children[[j]]) == 1L) return(invisible())
    nxt <- (1 - p[, j]) * acc
    for (k in tree$children[[j]]) walk(k, nxt)
  }
  walk(tree$n_node, 1)
  p2
}

## posterior mu/sigma per leaf, mixing node columns with complement columns
leaf_posteriors <- function(tree, mu_full, sigma_full, p2, p.t) {
  n <- nrow(mu_full); n_node <- tree$n_node
  tail_cols <- seq.int(n_node + 1L, tree$n_col)
  mu_post <- sigma_post <- matrix(0, n, tree$n_leaf)
  paths <- tree_paths(tree)
  for (i in seq_len(tree$n_leaf)) {
    path <- paths[[i]]
    index <- c(path, tail_cols)
    ## a complement of an ancestor carries the same information as the ancestor
    hit <- which(tree$comp_of %in% path)
    if (length(hit)) index[length(path) + hit] <- tree$comp_of[hit]
    weight <- cbind(p2[, path, drop = FALSE], p.t[, -1, drop = FALSE])
    mu_post[, i]    <- rowSums(mu_full[, index, drop = FALSE]    * weight)
    sigma_post[, i] <- rowSums(sigma_full[, index, drop = FALSE] * weight)
  }
  list(mu = mu_post, sigma = sigma_post)
}

## ---- one-vs-rest markers ------------------------------------------------------
## Vectorised over CpGs by grouping on the winning node.
##
## p.t column k+1 is complement position k, which refers to node comp_of[k];
## that node's own column is `node` and its complement column is n_node + k.
select_onevrest <- function(tree, m, p.t, mu_full, sigma_full, IQR, p_post, Nmean) {
  n <- nrow(p.t)
  n_comp <- length(tree$comp_of)
  strict_ct <- unique(c(as.integer(unlist(m$strict_celltypes)),
                        which(tree$leaf_sizes <= m$strict_max_samples)))
  k  <- max.col(p.t, ties.method = "first") - 1L
  mx <- p.t[cbind(seq_len(n), k + 1L)]
  cand <- which(!is.na(mx) & k >= 1L & k <= n_comp &
                mx >= m$lowest_p & Nmean >= m$min_read_depth)
  if (!length(cand)) return(NULL)

  out <- list()
  for (kk in sort(unique(k[cand]))) {
    node <- tree$comp_of[kk]
    tgt  <- tree$leaf[[node]]
    if (length(tgt) >= tree$n_leaf) next
    idx  <- cand[k[cand] == kk]
    ccol <- tree$n_node + kk

    rel <- abs((mu_full[idx, ccol] - mu_full[idx, node]) /
               sqrt(sigma_full[idx, ccol]^2 + sigma_full[idx, node]^2))
    thr <- if (node %in% strict_ct) m$strict_rel_distance else m$min_rel_distance

    keep <- !is.na(rel) & rel >= thr
    keep <- keep & rowSums(IQR[idx, tgt,  drop = FALSE] > m$min_IQR_target,    na.rm = TRUE) == 0
    keep <- keep & rowSums(IQR[idx, -tgt, drop = FALSE] > m$min_IQR_offtarget, na.rm = TRUE) == 0
    keep <- keep & abs(rowMeans(p_post[idx, tgt,  drop = FALSE]) -
                       rowMeans(p_post[idx, -tgt, drop = FALSE])) >= m$min_abs_distance
    if (!any(keep)) next
    out[[length(out) + 1L]] <- data.frame(
      idx = idx[keep],
      celltype = if (node <= tree$n_leaf) tree$celltypes[node] else as.character(node),
      node = node,
      ## evidence strength, used when a per-label cap has to choose
      score = mx[idx[keep]], stringsAsFactors = FALSE)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}

## ---- semi-specific markers ------------------------------------------------------
select_semi <- function(tree, m, p.t, p2, p_split, mu_post, sigma_post, IQR, Nmean, p_post) {
  n_leaf <- tree$n_leaf
  strict_ct <- as.integer(unlist(m$strict_celltypes))
  gate_thr  <- m$lowest_p
  pair_mode <- identical(m$comparison %||% "target_vs_rest", "pair")
  want_sib  <- isTRUE(m$gate_requires_sibling_block)
  sib_mode  <- m$sibling_in_comparison %||% "counted"
  rel_thr0  <- m$min_rel_distance
  use_rel   <- !is.null(rel_thr0) && !is.na(rel_thr0)
  eff       <- m$min_abs_distance_T0
  slack     <- m$slack

  sibling_of <- function(v) { pa <- tree$parent[v]
    if (pa == v) return(integer(0)); setdiff(tree$children[[pa]], v) }
  gate_for <- function(node) {
    g <- p2[, node]
    if (want_sib) { sb <- sibling_of(node); if (length(sb)) g <- g * p_split[, sb[1]] }
    g
  }
  base <- which(p.t[, 1] >= m$lowest_p & Nmean >= m$min_read_depth)
  out <- list()
  add <- function(rows, node) if (length(rows))
    out[[length(out) + 1L]] <<- data.frame(idx = rows, celltype = as.character(node),
                                           node = node, score = p2[rows, node],
                                           stringsAsFactors = FALSE)

  for (node in seq_len(tree$n_node - 1L)) {
    tgt <- tree$leaf[[node]]
    if (length(tgt) >= n_leaf) next
    g <- gate_for(node)
    idx <- base[!is.na(g[base]) & g[base] >= gate_thr]
    if (!length(idx)) next

    sb <- sibling_of(node)
    sib_leaves <- if (length(sb)) tree$leaf[[sb[1]]] else integer(0)
    cmp <- setdiff(seq_len(n_leaf), tgt)
    if (pair_mode && identical(sib_mode, "excluded")) cmp <- setdiff(cmp, sib_leaves)
    if (!length(cmp)) next
    is_strict <- node %in% strict_ct
    need <- if (!pair_mode) {
      if (is_strict) n_leaf - slack else n_leaf - length(tgt) - slack
    } else if (identical(sib_mode, "exempt")) {
      length(cmp) - length(sib_leaves) - slack
    } else length(cmp) - slack
    if (need < 1L) next
    rel_thr <- if (is_strict) m$strict_rel_distance else rel_thr0

    mu <- mu_post[idx, , drop = FALSE]; sg <- sigma_post[idx, , drop = FALSE]
    px <- p_post[idx, , drop = FALSE]
    mt <- rowMeans(mu[, tgt, drop = FALSE]); st <- rowMeans(sg[, tgt, drop = FALSE])
    pt <- rowMeans(px[, tgt, drop = FALSE])

    keep <- rep(TRUE, length(idx))
    if (use_rel) {
      dd <- abs((if (pair_mode) mu[, cmp, drop = FALSE] else mu) - mt) /
            sqrt(st^2 + (if (pair_mode) sg[, cmp, drop = FALSE] else sg)^2)
      keep <- keep & rowSums(dd > rel_thr, na.rm = TRUE) >= need
    }
    keep <- keep &
      rowSums(IQR[idx, tgt, drop = FALSE] >= m$min_IQR_target, na.rm = TRUE) == 0
    if (pair_mode && length(sib_leaves))
      keep <- keep &
        rowSums(IQR[idx, sib_leaves, drop = FALSE] > m$min_IQR_target, na.rm = TRUE) == 0
    bg <- if (pair_mode) setdiff(seq_len(n_leaf), c(tgt, sib_leaves)) else -tgt
    if (length(bg))
      keep <- keep & rowSums(IQR[idx, bg, drop = FALSE] < m$min_IQR_offtarget,
                             na.rm = TRUE) >= (if (pair_mode) length(bg) - slack else need + 2L)
    keep <- keep &
      rowSums(abs(px[, cmp, drop = FALSE] - pt) > eff, na.rm = TRUE) >=
        (if (pair_mode) need else need + 1L)
    if (any(keep)) add(idx[keep], node)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

## ---- sibling-contrast markers -----------------------------------------------
select_sibling_contrast <- function(tree, m, p, p.t, mu_post, sigma_post, IQR, Nmean, p_post) {
  pw <- m$sibling_contrast
  if (!isTRUE(pw$enabled)) return(NULL)
  bg_mode  <- pw$background %||% "iqr"
  if (!bg_mode %in% c("iqr", "joint_p", "span", "both", "none"))
    die("markers.sibling_contrast.background must be one of iqr, joint_p, span, both, none")
  max_span <- pw$max_bg_span %||% 0.25
  n_leaf <- tree$n_leaf
  n      <- nrow(mu_post)
  slack  <- m$slack

  anc_of <- function(x) { o <- integer(0)
    while (x != tree$n_node) { x <- tree$parent[x]; o <- c(o, x) }; o }
  pair_joint <- function(j) {
    g <- tree$children[[j]]; acc <- rep(1, n)
    for (k in anc_of(j)) acc <- acc * (1 - p[, k])
    acc * (1 - p[, j]) * p[, g[1]] * p[, g[2]] * p.t[, 1]
  }
  rowMean <- function(M, j)
    if (length(j) == 1L) M[, j] else rowMeans(M[, j, drop = FALSE])

  depth_ok <- Nmean >= m$min_read_depth
  out <- list()
  for (j in which(lengths(tree$children) == 2L)) {
    g <- tree$children[[j]]
    A <- tree$leaf[[g[1]]]; B <- tree$leaf[[g[2]]]
    R <- setdiff(seq_len(n_leaf), c(A, B))
    if (!length(R)) next
    mA <- rowMean(mu_post, A);    mB <- rowMean(mu_post, B)
    sA <- rowMean(sigma_post, A); sB <- rowMean(sigma_post, B)
    pA <- rowMean(p_post, A);     pB <- rowMean(p_post, B)
    rel <- abs(mA - mB) / sqrt(sA^2 + sB^2)
    keep <- depth_ok & !is.na(rel) &
      rel >= m$min_rel_distance &
      abs(pA - pB) >= m$min_abs_distance &
      rowSums(IQR[, A, drop = FALSE] >= m$min_IQR_target, na.rm = TRUE) == 0 &
      rowSums(IQR[, B, drop = FALSE] >= m$min_IQR_target, na.rm = TRUE) == 0
    if (bg_mode %in% c("joint_p", "both"))
      keep <- keep & pair_joint(j) >= m$lowest_p
    if (bg_mode %in% c("span", "both")) {
      span <- if (length(R) > 1L)
        apply(p_post[, R, drop = FALSE], 1, function(v) diff(range(v))) else rep(0, n)
      keep <- keep & span <= max_span
    }
    if (identical(bg_mode, "iqr"))
      keep <- keep &
        rowSums(IQR[, R, drop = FALSE] < m$min_IQR_offtarget, na.rm = TRUE) >=
          max(0L, length(R) - slack)
    rows <- which(keep)
    if (!length(rows)) next
    sc <- pair_joint(j)[rows]
    for (side in g)
      out[[length(out) + 1L]] <- data.frame(idx = rows, celltype = as.character(side),
                                            node = side, score = sc,
                                            stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

## Alias retained for callers that predate the rename.
select_T0 <- function(tree, m, p.t, p2, mu_post, sigma_post, IQR, Nmean, p_post)
  select_semi(tree, m, p.t, p2, p2, mu_post, sigma_post, IQR, Nmean, p_post)

## ---- driver -------------------------------------------------------------------
select_markers <- function(tree, cfg_m, ll_full, ll_trunc, mu_full, sigma_full,
                           pi.j, pi.t, Nmean, cpg_index) {
  ## Fail before the expensive part rather than silently ignoring a stale key.
  if (!is.null(cfg_m$pairwise))
    die("markers.pairwise was renamed markers.sibling_contrast -- rename it in ",
        "the config rather than leaving a key that is silently ignored")
  tree$leaf_sizes <- lengths(tree$leaf)[seq_len(tree$n_leaf)]
  pa  <- calculate.p.gj(tree$children, pi.j, ll_full, tree$layer)
  p.t <- calculate.p.t(pi.t, pi.j, tree, ll_full, ll_trunc, pa[[2]])
  p2  <- node_joint_p(tree, pa[[1]]) * p.t[, 1]

  lp  <- leaf_posteriors(tree, mu_full, sigma_full, p2, p.t)
  IQR <- beta_IQR(lp$mu, lp$sigma, cfg_m$alpha_cap)
  p_post <- pnorm(lp$mu / sqrt(1 + lp$sigma^2))
  ## Cross-cell-type median IQR: how variable this CpG is within a cell type,
  ## summarised over cell types.  The per-label cap ranks on it, smallest first.
  median_iqr <- apply(IQR, 1, median, na.rm = TRUE)

  log_msg("  one-vs-rest scan (celltype and class) ...")
  ovr <- select_onevrest(tree, cfg_m, p.t, mu_full, sigma_full, IQR, p_post, Nmean)
  log_msg("  semi-specific scan ...")
  t0  <- select_semi(tree, cfg_m, p.t, p2, pa[[1]], lp$mu, lp$sigma, IQR, Nmean, p_post)

  pw <- NULL
  if (isTRUE(cfg_m$sibling_contrast$enabled)) {
    log_msg("  sibling-contrast scan ...")
    pw <- select_sibling_contrast(tree, cfg_m, pa[[1]], p.t, lp$mu, lp$sigma, IQR, Nmean, p_post)
  }

  if (!is.null(ovr)) ovr$onevsrest <- 1L
  if (!is.null(t0))  t0$onevsrest  <- 0L
  if (!is.null(pw))  pw$onevsrest  <- 0L
  sig <- rbind(ovr, t0, pw)
  if (is.null(sig) || !nrow(sig)) { log_msg("  no markers passed"); return(NULL) }
  ## One row per (CpG, node, kind).  The two scans above cannot collide -- each
  ## visits a node once -- so this is a no-op unless the sibling-contrast scan ran and
  ## re-found a node that select_semi had already reported for the same CpG.
  sig <- sig[!duplicated(sig[, c("idx", "node", "onevsrest")]), , drop = FALSE]
  sig <- sig[order(sig$idx), , drop = FALSE]

  sig$semi_pair <- NA
  s <- which(sig$onevsrest == 0L)
  if (length(s)) {
    sib <- vapply(sig$node[s], function(v)
      as.integer(setdiff(tree$children[[tree$parent[v]]], v)[1]), 0L)
    sig$semi_pair[s] <- p2[cbind(sig$idx[s], sig$node[s])] *
                        pa[[1]][cbind(sig$idx[s], sib)] >= cfg_m$lowest_p
  }

  mu_cols <- lp$mu[sig$idx, , drop = FALSE]
  sd_cols <- lp$sigma[sig$idx, , drop = FALSE]
  colnames(mu_cols) <- paste0("mu_",    tree$celltypes)
  colnames(sd_cols) <- paste0("sigma_", tree$celltypes)
  out <- cbind(data.frame(index = cpg_index[sig$idx], local_idx = sig$idx,
                          celltype = sig$celltype, node = sig$node,
                          onevsrest = sig$onevsrest, semi_pair = sig$semi_pair,
                          score = sig$score,
                          median_iqr = median_iqr[sig$idx],
                          stringsAsFactors = FALSE),
               mu_cols, sd_cols)
  nmulti <- sum(table(out$index) > 1L)
  is_leaf <- out$node <= tree$n_leaf
  log_msg("  markers: ", nrow(out), " rows over ", length(unique(out$index)),
          " CpGs (celltype ", sum(out$onevsrest == 1L & is_leaf),
          ", class ", sum(out$onevsrest == 1L & !is_leaf),
          ", semi_pair ", sum(out$onevsrest == 0L & out$semi_pair %in% TRUE),
          ", semi_single ", sum(out$onevsrest == 0L & out$semi_pair %in% FALSE),
          "); ", nmulti,
          " CpGs resolve more than one block")
  out
}

## cell-type-specific beta matrix from the posterior mu/sigma columns
signature_beta <- function(sig, celltypes) {
  mu <- as.matrix(sig[, paste0("mu_",    celltypes), drop = FALSE])
  sg <- as.matrix(sig[, paste0("sigma_", celltypes), drop = FALSE])
  b  <- pnorm(mu / sqrt(1 + sg^2))
  colnames(b) <- celltypes
  b
}
