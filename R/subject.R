## cf-TREBLE pipeline -- subject-specific signature ----------------------------
##
## Replaces each cell type's single value with one per reference subject:
##
##   p_subj[g, i] = sum over partitions t  P(t | g) * p_Ind[g, side_of(t, i), i]
##
## A partition is either T0 (the tree itself, weight p2 per node) or Tk (node k
## against its complement, weight p.t[k+1]).  For sample i only ONE side of each
## partition contains it, and that is the side whose PostExp_IndSig value counts.
## So the weight attached to node n, for a sample in leaf j, is
##
##   n on j's path to the root : p2[n] + p.t[complement of n]   (i is on n's side)
##   n = complement of some b  : p.t[n]                          (i is in n iff b
##                                                                is not on the path)
##
## The weights sum to one, and every p_Ind value is in [0, 1]; dropping a set of
## nodes whose weights sum to w therefore moves the answer by at most w.  That
## is what `weight_coverage` trades away, and it is the whole speedup.

## For each leaf j, the node-indexed weight map.  Returns a list with, per leaf,
## the node ids and, for each, which of the two weight sources feed it.
subject_weight_plan <- function(tree) {
  paths <- tree_paths(tree)
  comp_col <- integer(tree$n_node)                 # node -> its complement column
  comp_col[tree$comp_of] <- tree$n_node + seq_along(tree$comp_of)
  lapply(seq_len(tree$n_leaf), function(j) {
    pth <- paths[[j]]
    ## complements that contain leaf j: those whose base node is off the path
    off <- which(!(tree$comp_of %in% pth))
    list(path = pth,
         path_comp = comp_col[pth],               # 0 where the node has no complement
         comp_nodes = tree$n_node + off,          # complement columns containing j
         comp_w = off + 1L)                       # their column in p.t
  })
}

## Weight vector over nodes for one CpG and one leaf.
subject_weights_one <- function(plan_j, p2_g, p.t_g) {
  n_node <- length(p2_g)
  wp <- p2_g[plan_j$path]
  hc <- plan_j$path_comp > 0L
  ## p.t column 1 is T0; complement column n_node+k is p.t column k+1
  if (any(hc)) wp[hc] <- wp[hc] + p.t_g[plan_j$path_comp[hc] - n_node + 1L]
  c(setNames(wp, plan_j$path), setNames(p.t_g[plan_j$comp_w], plan_j$comp_nodes))
}

## Which nodes must be evaluated for this CpG, given a weight-coverage target.
## Union over leaves, because PostExp_IndSig runs per node over all its samples.
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

## p_Ind for one CpG over a given set of nodes.
## Returns a list: for each node, a named vector of per-sample posteriors.
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

## One CpG -> one row of p_subj (length = number of reference samples).
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
  ## Renormalise by the weight actually used.  Leaving it out would implicitly
  ## score every dropped partition as p = 0; dividing instead treats them as the
  ## weighted mean of the kept ones, which is the better guess and shrinks the
  ## error from (1 - coverage) to (1 - coverage) x |p_dropped - p_kept|.
  ## With coverage = 1 the divisor is 1, so the exact path is untouched.
  ## A sample with no reads at this CpG has no weight at all -> NA, not 0.
  ok <- wsum > 1e-12
  out[ok] <- out[ok] / wsum[ok]
  out[!ok] <- NA_real_
  list(p = out, covered = wsum, n_nodes = length(nodes))
}
