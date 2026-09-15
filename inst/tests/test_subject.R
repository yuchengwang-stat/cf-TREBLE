#!/usr/bin/env Rscript
## Checks on the subject-specific weighting, no pipeline run needed.
##   1. the weights for every sample sum to 1 (they are a convex combination)
##   2. they agree with an independent construction of the same weights
##   3. pruning to coverage c really leaves at least c of the mass, per sample
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "..", "scripts", "boot.R"))
a <- args_parse(list(tree_lists = "", reps = "40", seed = "1"))
set.seed(as.integer(a$seed)); reps <- as.integer(a$reps)

if (nzchar(a$tree_lists)) {
  tree <- tree_from_lists(file.path(a$tree_lists, "binary.tree2.0.RDS"),
                          file.path(a$tree_lists, "ancestor.tree2.0.RDS"),
                          file.path(a$tree_lists, "all_celltype_sample_2.0.RDS"))
} else {
  nl <- 16L; hc <- hclust(dist(matrix(rnorm(nl * 30), nl)), "ward.D2")
  ctn <- sprintf("CT%02d", seq_len(nl)); hc$labels <- ctn
  tree <- hclust_to_tree(hc, ctn)
}
plan  <- subject_weight_plan(tree)
paths <- tree_paths(tree)
cat(sprintf("tree: %d leaves, %d nodes, %d columns\n", tree$n_leaf, tree$n_node, tree$n_col))

## A draw of (p2, p.t) with the structure the real quantities have.  p2 must be
## built the way node_joint_p builds it -- one unit of mass walking down the
## tree -- so that along ANY root-to-leaf path the p2 values sum to exactly
## p.t[1].  (Leaves have split probability 1, which is what closes the sum.)
draw <- function() {
  p.t <- c(rgamma(1, 9), rgamma(length(tree$comp_of), 0.05)); p.t <- p.t / sum(p.t)
  pmat <- matrix(runif(tree$n_node, 0.05, 0.95), 1, tree$n_node)
  pmat[1, seq_len(tree$n_leaf)] <- 1                       # a leaf cannot split further
  p2 <- drop(node_joint_p(tree, pmat)) * p.t[1]
  list(p2 = p2, p.t = p.t)
}
## the same weights built by index arithmetic instead of the plan
index_form_w <- function(j, p2, p.t) {
  pth <- paths[[j]]
  index <- c(pth, seq.int(tree$n_node + 1L, tree$n_col))
  hit <- which(tree$comp_of %in% pth)
  index[length(pth) + hit] <- tree$comp_of[hit]
  tapply(c(p2[pth], p.t[-1]), index, sum)
}

pass <- TRUE
A <- function(nm, ok, d = "") { cat(sprintf("  %-44s %s %s\n", nm, if (isTRUE(ok)) "PASS" else "FAIL", d))
                                pass <<- pass && isTRUE(ok) }

cat("\n=== 1. weights are a convex combination ===\n")
worst <- 0
for (r in seq_len(reps)) { d <- draw()
  for (j in seq_len(tree$n_leaf))
    worst <- max(worst, abs(sum(subject_weights_one(plan[[j]], d$p2, d$p.t)) - 1)) }
A("weights sum to 1 for every leaf", worst < 1e-9, sprintf("(max |sum-1| %.2e)", worst))

cat("\n=== 2. identical to the index-arithmetic construction ===\n")
maxd <- 0
for (r in seq_len(min(reps, 20))) { d <- draw()
  for (j in seq_len(tree$n_leaf)) {
    mine <- subject_weights_one(plan[[j]], d$p2, d$p.t); leg <- index_form_w(j, d$p2, d$p.t)
    if (!setequal(names(mine), names(leg))) { maxd <- Inf; break }
    maxd <- max(maxd, max(abs(mine[names(leg)] - leg)))
  } }
A("same node set and same weights", maxd < 1e-12, sprintf("(max diff %.2e)", maxd))

cat("\n=== 3. pruning delivers the coverage it promises ===\n")
for (cov in c(0.99, 0.999, 0.9999)) {
  worst_cov <- 1; kept <- c()
  for (r in seq_len(min(reps, 25))) { d <- draw()
    nodes <- subject_nodes_needed(plan, d$p2, d$p.t, cov); kept <- c(kept, length(nodes))
    for (j in seq_len(tree$n_leaf)) {
      w <- subject_weights_one(plan[[j]], d$p2, d$p.t)
      worst_cov <- min(worst_cov, sum(w[as.integer(names(w)) %in% nodes]) / sum(w)) } }
  A(sprintf("coverage >= %.4f for every sample", cov), worst_cov >= cov - 1e-12,
    sprintf("(worst %.6f, nodes kept %d/%d)", worst_cov, round(median(kept)), tree$n_col))
}
cat("\n", if (pass) "SUBJECT TESTS PASSED" else "SUBJECT TESTS FAILED", "\n")
quit(status = if (pass) 0L else 1L)
