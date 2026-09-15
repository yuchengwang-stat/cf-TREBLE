## cf-TREBLE pipeline -- tree construction, validation, serialisation ----------
##
## Everything downstream reads ONE object, produced here.  No stage is allowed
## to hardcode a node count, a layer range, or a complement offset again.
##
## tree bundle fields
##   celltypes  chr[n_leaf]   leaf names, in column order of M/N
##   children   list[n_node]  internal -> c(left, right); leaf -> its own index
##   parent     int[n_node]   parent of each node; root -> itself
##   leaf       list[n_node]  set of leaves under each node
##   leaf_all   list[n_col]   leaf, then the complement of every node whose
##                            complement is not already some node's leaf set
##   comp_of    int[n_comp]   leaf_all[[n_node+k]] is the complement of node comp_of[k]
##   layer      list          internal nodes grouped by height, bottom-up
##   n_leaf / n_node / n_col

set_key <- function(v) paste(sort(as.integer(v)), collapse = ",")

## leaf sets, bottom-up from `children`
tree_leaf_sets <- function(children) {
  n <- length(children)
  leaf <- vector("list", n)
  repeat {
    done <- TRUE
    for (i in seq_len(n)) {
      if (!is.null(leaf[[i]])) next
      kid <- children[[i]]
      if (length(kid) == 1L) { leaf[[i]] <- as.integer(i); done <- FALSE; next }
      if (all(!vapply(leaf[kid], is.null, TRUE))) {
        leaf[[i]] <- sort(unlist(leaf[kid], use.names = FALSE)); done <- FALSE
      }
    }
    if (done) break
  }
  if (any(vapply(leaf, is.null, TRUE))) die("children list is not a rooted tree")
  leaf
}

## height above the leaves; leaves are 0
tree_heights <- function(children) {
  n <- length(children); h <- rep(NA_integer_, n)
  repeat {
    done <- TRUE
    for (i in seq_len(n)) {
      if (!is.na(h[i])) next
      kid <- children[[i]]
      if (length(kid) == 1L) { h[i] <- 0L; done <- FALSE; next }
      if (all(!is.na(h[kid]))) { h[i] <- max(h[kid]) + 1L; done <- FALSE }
    }
    if (done) break
  }
  h
}

## `layer`: internal nodes grouped by height, bottom-up.  A node is only
## evaluated after both children, which is exactly what the EM recursion needs.
tree_layers <- function(children) {
  h <- tree_heights(children)
  internal <- which(lengths(children) == 2L)
  unname(split(internal, h[internal]))
}

tree_parent_from_children <- function(children) {
  n <- length(children); p <- seq_len(n)
  for (i in seq_len(n)) { kid <- children[[i]]; if (length(kid) == 2L) p[kid] <- i }
  p
}

## Complements, de-duplicated.  The complement of the root is empty, and the
## root's two children are each other's complement, so their one-vs-rest
## partition is the tree's own split at the root and carries nothing new.
tree_complements <- function(leaf) {
  n <- length(leaf)
  universe <- leaf[[n]]
  seen <- vapply(leaf, set_key, "")
  comp <- list(); comp_of <- integer(0)
  for (i in seq_len(n)) {
    cc <- setdiff(universe, leaf[[i]])
    if (!length(cc)) next
    k <- set_key(cc)
    if (k %in% seen) next
    seen <- c(seen, k)
    comp[[length(comp) + 1L]] <- cc
    comp_of <- c(comp_of, i)
  }
  list(comp = comp, comp_of = comp_of)
}

build_tree <- function(children, celltypes, parent = NULL) {
  children <- lapply(children, as.integer)
  n_node <- length(children)
  leaf <- tree_leaf_sets(children)
  n_leaf <- length(leaf[[n_node]])
  if (length(celltypes) != n_leaf)
    die("celltypes has ", length(celltypes), " entries but the tree has ", n_leaf, " leaves")
  cp <- tree_complements(leaf)
  if (is.null(parent)) parent <- tree_parent_from_children(children)
  parent <- as.integer(unlist(parent))
  tr <- list(celltypes = as.character(celltypes),
             children = children, parent = parent, leaf = leaf,
             leaf_all = c(leaf, cp$comp), comp_of = cp$comp_of,
             layer = tree_layers(children),
             n_leaf = n_leaf, n_node = n_node, n_col = n_node + length(cp$comp))
  validate_tree(tr)
  tr
}

validate_tree <- function(tr) {
  n <- tr$n_node
  if (tr$n_node != 2L * tr$n_leaf - 1L)
    die("not a binary tree: n_node=", tr$n_node, " n_leaf=", tr$n_leaf)
  if (!identical(sort(which(lengths(tr$children) == 1L)), seq_len(tr$n_leaf)))
    die("leaves must be nodes 1..n_leaf, in order")
  if (tr$parent[n] != n) die("root (node ", n, ") must be its own parent")
  ## every non-root node reaches the root
  for (i in seq_len(n - 1L)) {
    cur <- i; steps <- 0L
    while (cur != n) { cur <- tr$parent[cur]; steps <- steps + 1L
      if (steps > n) die("parent chain from node ", i, " does not terminate") }
  }
  ## children/parent agree
  for (i in seq_len(n)) for (k in tr$children[[i]])
    if (length(tr$children[[i]]) == 2L && tr$parent[k] != i)
      die("parent/children disagree at node ", i)
  ## layers cover every internal node exactly once, in a valid order
  seen <- unlist(tr$layer)
  if (!setequal(seen, which(lengths(tr$children) == 2L)) || anyDuplicated(seen))
    die("layer does not partition the internal nodes")
  invisible(TRUE)
}

## path from a leaf up to the root (inclusive), used to weight node posteriors
tree_paths <- function(tr) {
  lapply(seq_len(tr$n_node), function(i) {
    path <- i
    while (path[length(path)] != tr$n_node) path <- c(path, tr$parent[path[length(path)]])
    path
  })
}

## ---- entry points ------------------------------------------------------------
tree_from_lists <- function(children_rds, parent_rds, celltypes_rds) {
  children <- readRDS(children_rds)
  parent   <- if (!is.null(parent_rds) && nzchar(parent_rds)) readRDS(parent_rds) else NULL
  ct       <- readRDS(celltypes_rds)
  build_tree(children, if (is.list(ct)) names(ct) else as.character(ct), parent)
}

## bed: rows = loci, cols = c(chr, start, end, samples...)
tree_from_bed <- function(bed, celltype_of_sample, knn_k = 5,
                          dist_method = "euclidean", linkage = "ward.D2") {
  if (!requireNamespace("impute", quietly = TRUE))
    die("package 'impute' (Bioconductor) is required for --from-bed")
  meta  <- bed[, 1:3, drop = FALSE]
  score <- as.matrix(bed[, -(1:3), drop = FALSE])
  score <- score[rowMeans(is.na(score)) <= 0.5, , drop = FALSE]
  imp <- impute::impute.knn(t(score), k = knn_k)$data       # features x samples
  hc  <- hclust(dist(imp, method = dist_method), method = linkage)
  list(hclust = hc, celltype_of_sample = celltype_of_sample)
}

## Build a cell-type tree straight from the reference counts.  Self-contained:
## no external bed, and the sample names are the ones the pipeline already uses,
## so it cannot disagree with sample_list.
##
##   m_chunks/n_chunks : sprintf patterns for the count matrices
##   chunks            : which to draw CpGs from
##   n_cpg             : how many CpGs to sample per chunk
tree_from_reference <- function(m_chunks, n_chunks, chunks, sample_list,
                                n_cpg = 20000, min_depth = 10, min_sd = 0.05,
                                dist_method = "euclidean", linkage = "ward.D2",
                                seed = 1) {
  groups <- rep(seq_along(sample_list), times = lengths(sample_list))
  prof <- list()
  for (k in chunks) {
    M <- readRDS(sprintf(m_chunks, k)); N <- readRDS(sprintf(n_chunks, k))
    set.seed(seed + k)
    rows <- sort(sample.int(nrow(M), min(n_cpg, nrow(M))))
    M <- M[rows, , drop = FALSE]; N <- N[rows, , drop = FALSE]
    Ms <- t(rowsum(t(M), group = groups)); Ns <- t(rowsum(t(N), group = groups))
    beta <- Ms / Ns
    ## every cell type must be measured, and the CpG must actually vary
    ok <- stats::complete.cases(beta) &
          apply(Ns, 1, min) >= min_depth &
          apply(beta, 1, sd) >= min_sd
    if (any(ok)) prof[[length(prof) + 1L]] <- beta[ok, , drop = FALSE]
    rm(M, N, Ms, Ns, beta); gc()
  }
  if (!length(prof)) die("no CpG passed the tree-building filters")
  B <- do.call(rbind, prof)
  colnames(B) <- names(sample_list)
  log_msg("tree from reference: ", nrow(B), " informative CpGs x ", ncol(B), " cell types")
  hc <- hclust(dist(t(B), method = dist_method), method = linkage)
  list(hclust = hc, n_cpg_used = nrow(B), profile = B)
}

## hclust over CELL TYPES -> children/parent in the 1..n_leaf, n_leaf+1..2n-1
## numbering the rest of the pipeline assumes.
hclust_to_tree <- function(hc, celltypes = hc$labels) {
  n_leaf <- length(hc$order)
  n_node <- 2L * n_leaf - 1L
  children <- vector("list", n_node)
  for (i in seq_len(n_leaf)) children[[i]] <- i
  for (m in seq_len(n_leaf - 1L)) {
    node <- n_leaf + m
    children[[node]] <- as.integer(vapply(hc$merge[m, ], function(x)
      if (x < 0) -x else n_leaf + x, 0L))
  }
  if (is.null(celltypes)) die("hclust has no labels; pass celltypes= explicitly")
  build_tree(children, celltypes)
}
