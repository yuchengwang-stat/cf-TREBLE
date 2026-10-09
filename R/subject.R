subject_weight_plan <- function(tree) {
  paths <- tree_paths(tree)
  comp_col <- integer(tree$n_node)
  comp_col[tree$comp_of] <- tree$n_node + seq_along(tree$comp_of)
  lapply(seq_len(tree$n_leaf), function(j) {
    pth <- paths[[j]]
    off <- which(!(tree$comp_of %in% pth))
    list(path = pth,
         path_comp = comp_col[pth],
         comp_nodes = tree$n_node + off,
         comp_w = off + 1L)
  })
}

subject_weights_one <- function(plan_j, p2_g, p.t_g) {
  n_node <- length(p2_g)
  wp <- p2_g[plan_j$path]
  hc <- plan_j$path_comp > 0L
  if (any(hc)) wp[hc] <- wp[hc] + p.t_g[plan_j$path_comp[hc] - n_node + 1L]
  c(setNames(wp, plan_j$path), setNames(p.t_g[plan_j$comp_w], plan_j$comp_nodes))
}

subject_nodes_needed <- function(plan, p2_g, p.t_g, coverage = 0.9999) {
  keep <- integer(0)
  for (pj in plan) {
    w <- subject_weights_one(pj, p2_g, p.t_g)
    o <- order(w, decreasing = TRUE)
    cum <- cumsum(w[o]) / sum(w)
    k <- which(cum >= coverage)[1]
    if (is.na(k)) k <- length(w)
    keep <- union(keep, as.integer(names(w)[o[seq_len(k)]]))
  }
  sort(keep)
}

subject_pInd_one <- function(nodes, tree, Mrow, Urow, groups, pr) {
  lapply(nodes, function(nd) {
    cols <- which(groups %in% tree$leaf_all[[nd]])
    v <- tryCatch(
      PostExp_IndSig(M = Mrow[cols], U = Urow[cols],
                     sigma.grid = pr$sigma.grid, sigma.pi = pr$sigma.pi,
                     mu.mu = pr$delta, mu.tau2 = pr$tau2, mu.pi = pr$lambda,
                     Gauss.grid = pr$gh.x, Gauss.w = pr$gh.w),
      error = function(e) rep(NA_real_, length(cols)))
    setNames(v, cols)
  })
}

subject_one_cpg <- function(tree, plan, groups, Mrow, Urow, p2_g, p.t_g, pr,
                            coverage = 0.9999) {
  nodes <- if (coverage >= 1) seq_len(tree$n_col)
           else subject_nodes_needed(plan, p2_g, p.t_g, coverage)
  pind <- subject_pInd_one(nodes, tree, Mrow, Urow, groups, pr)
  names(pind) <- as.character(nodes)

  out <- numeric(length(groups)); wsum <- numeric(length(groups))
  for (i in seq_along(groups)) {
    w <- subject_weights_one(plan[[groups[i]]], p2_g, p.t_g)
    w <- w[names(w) %in% names(pind)]
    if (!length(w)) next
    for (nm in names(w)) {
      v <- pind[[nm]][[as.character(i)]]
      if (!is.null(v) && !is.na(v)) { out[i] <- out[i] + w[[nm]] * v; wsum[i] <- wsum[i] + w[[nm]] }
    }
  }
  ok <- wsum > 1e-12
  out[ok] <- out[ok] / wsum[ok]
  out[!ok] <- NA_real_
  list(p = out, covered = wsum, n_nodes = length(nodes))
}
