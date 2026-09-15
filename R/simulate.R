## cf-TREBLE pipeline -- synthetic reference and plasma, with known truth ------
##
## Used by bin/validate.R to answer "are the numbers right", as opposed to
## bin/selftest.R which only answers "does the wiring hold".

## Cell-type-level methylation truth, built ON the tree so that both kinds of
## marker the pipeline looks for actually exist:
##   * one-vs-rest markers: one leaf differs from all the others
##   * class markers:       an internal node's whole leaf set differs
simulate_truth <- function(tree, n_cpg, marker_frac = 0.30,
                           class_frac = 0.5, sep = 0.80) {
  n_leaf <- tree$n_leaf
  P <- matrix(runif(n_cpg * n_leaf, 0.25, 0.75), n_cpg, n_leaf)   # background
  kind <- rep("bg", n_cpg); node_of <- rep(NA_integer_, n_cpg)

  n_mark <- round(n_cpg * marker_frac)
  if (n_mark > 0) {
    rows <- sample.int(n_cpg, n_mark)
    n_cls <- round(n_mark * class_frac)
    internal <- which(lengths(tree$children) == 2L)
    internal <- internal[vapply(internal, function(j)
      length(tree$leaf[[j]]) < n_leaf - 1L, TRUE)]           # must leave an "off" side
    for (i in seq_along(rows)) {
      g <- rows[i]
      if (i <= n_cls && length(internal)) {
        node <- sample(internal, 1); kind[g] <- "class"
      } else {
        node <- sample.int(n_leaf, 1); kind[g] <- "onevrest"
      }
      node_of[g] <- node
      tgt <- tree$leaf[[node]]
      hi <- runif(1) < 0.5
      P[g, ]    <- if (hi) runif(n_leaf, 0, 1 - sep) else runif(n_leaf, sep, 1)
      P[g, tgt] <- if (hi) runif(length(tgt), sep, 1) else runif(length(tgt), 0, 1 - sep)
    }
  }
  list(P = P, kind = kind, node_of = node_of)
}

## Reference counts.  Samples of a cell type scatter around its truth, which is
## what the subject-level sigma in the model is there to absorb.
simulate_reference <- function(truth, groups, depth = 40, sample_sd = 0.05) {
  n_cpg <- nrow(truth$P); nS <- length(groups)
  Ps <- truth$P[, groups, drop = FALSE]
  Ps <- pmin(pmax(Ps + matrix(rnorm(n_cpg * nS, 0, sample_sd), n_cpg, nS), 1e-4), 1 - 1e-4)
  N <- matrix(rpois(n_cpg * nS, depth) + 5L, n_cpg, nS)
  list(M = matrix(rbinom(n_cpg * nS, N, Ps), n_cpg, nS), N = N, P_sample = Ps)
}

## A plasma sample at known cell-type proportions, written as a .beta file.
simulate_plasma <- function(truth, props, depth, path) {
  stopifnot(abs(sum(props) - 1) < 1e-8, length(props) == ncol(truth$P))
  n_cpg <- nrow(truth$P)
  p_mix <- as.vector(truth$P %*% props)
  N <- rpois(n_cpg, depth) + 1L
  M <- rbinom(n_cpg, N, p_mix)
  writeBin(as.integer(t(cbind(pmin(M, 255L), pmin(N, 255L)))), path, size = 1L)
  path
}

## ---- mixtures worth testing ----------------------------------------------------
make_mixtures <- function(celltypes) {
  n <- length(celltypes)
  z <- function() setNames(numeric(n), celltypes)
  mk <- function(v) { v <- pmax(v, 0); v / sum(v) }
  blood <- grep("^Blood|Megakary|Eryth", celltypes)
  if (!length(blood)) blood <- seq_len(min(5, n))
  other <- setdiff(seq_len(n), blood)

  ## 1. blood-dominated, like real cfDNA
  m1 <- z(); m1[blood] <- c(0.55, 0.20, rep(0.05, max(0, length(blood) - 2)))[seq_along(blood)]
  m1[sample(other, 3)] <- c(0.06, 0.03, 0.01)
  ## 2/2b. a small spike on a blood background -- the detection-limit case
  mk_spike <- function(frac) {
    v <- z(); v[blood] <- 1 / length(blood); v <- mk(v) * (1 - frac)
    v[spike] <- frac; structure(v, spike = celltypes[spike], spike_frac = frac)
  }
  spike <- sample(other, 1)
  m2  <- mk_spike(0.01)
  m2b <- mk_spike(0.001)
  ## 3. flat over ten random types
  m3 <- z(); m3[sample(seq_len(n), 10)] <- 0.1
  ## 4. two types only
  m4 <- z(); m4[sample(seq_len(n), 2)] <- c(0.7, 0.3)

  list("blood-dominated" = mk(m1), "1% spike" = m2, "0.1% spike" = m2b,
       "flat over 10" = mk(m3), "two types" = mk(m4))
}

## ---- accuracy ------------------------------------------------------------------
deconv_metrics <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  present <- truth > 0
  c(RMSE = sqrt(mean((est - truth)^2)),
    MAE = mean(abs(est - truth)),
    max_abs_err = max(abs(est - truth)),
    pearson = suppressWarnings(cor(est, truth)),
    spearman = suppressWarnings(cor(est, truth, method = "spearman")),
    err_on_present = if (any(present)) mean(abs(est[present] - truth[present])) else NA_real_,
    false_pos_mass = sum(est[!present]))
}
