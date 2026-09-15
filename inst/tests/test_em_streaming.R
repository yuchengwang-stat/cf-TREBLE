#!/usr/bin/env Rscript
## The streaming EM must give the same answer as the batch EM, not a close one.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "..", "scripts", "boot.R"))
a <- args_parse(list(n = "1200", parts = "5", seed = "3", tree_lists = ""))
set.seed(as.integer(a$seed))
if (nzchar(a$tree_lists)) {
  tree <- tree_from_lists(file.path(a$tree_lists, "binary.tree2.0.RDS"),
                          file.path(a$tree_lists, "ancestor.tree2.0.RDS"),
                          file.path(a$tree_lists, "all_celltype_sample_2.0.RDS"))
} else {
  nl <- 12L; hc <- hclust(dist(matrix(rnorm(nl * 30), nl)), "ward.D2")
  ctn <- sprintf("CT%02d", seq_len(nl)); hc$labels <- ctn; tree <- hclust_to_tree(hc, ctn)
}
n <- as.integer(a$n); parts <- as.integer(a$parts)
## Log-likelihoods must scale like real ones: a node's value is the likelihood
## of the samples under it, so it is on the order of the SUM of its children's.
## Independent columns make a parent ~half its children's total, which drives p
## to 1, collapses the weights, and makes the EM degenerate for reasons that
## have nothing to do with how the sum is accumulated.
leaf_ll <- matrix(rnorm(n * tree$n_leaf, -3, 0.6), n, tree$n_leaf)
build <- function(sets) vapply(sets, function(v)
  rowSums(leaf_ll[, v, drop = FALSE]), numeric(n))
ll  <- build(tree$leaf_all) + matrix(rnorm(n * tree$n_col, 0, 0.3), n, tree$n_col)
llt <- ll - abs(matrix(rnorm(n * tree$n_col, 0, 0.15), n, tree$n_col))

fit_b <- estimate_tree_prior(tree, ll, llt, tol = 1e-8, max_iter = 120)
grp <- split(seq_len(n), cut(seq_len(n), parts, labels = FALSE))
fit_s <- estimate_tree_prior_streaming(tree, seq_len(parts),
  function(k) list(ll = ll[grp[[k]], , drop = FALSE], llt = llt[grp[[k]], , drop = FALSE]),
  tol = 1e-8, max_iter = 120, progress = FALSE)

dj <- max(abs(fit_b$pi.j - fit_s$pi.j)); dt <- max(abs(fit_b$pi.t - fit_s$pi.t))
cat(sprintf("tree %d leaves; %d CpGs split into %d chunks\n", tree$n_leaf, n, parts))
cat(sprintf("  batch     : %d iters, final err %.3e\n", length(fit_b$trace), tail(fit_b$trace,1)))
cat(sprintf("  streaming : %d iters, final err %.3e\n", length(fit_s$trace), tail(fit_s$trace,1)))
cat(sprintf("  max |pi.j diff| %.3e\n  max |pi.t diff| %.3e\n", dj, dt))
ok <- dj < 1e-12 && dt < 1e-12 && length(fit_b$trace) == length(fit_s$trace)
cat("\n", if (ok) "STREAMING EM MATCHES BATCH" else "MISMATCH", "\n")
quit(status = if (ok) 0L else 1L)
