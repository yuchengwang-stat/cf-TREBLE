#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "..", "scripts", "boot.R"))
set.seed(5)

nl <- 12L
hc <- hclust(dist(matrix(rnorm(nl * 30), nl)), "ward.D2")
ctn <- sprintf("CT%02d", seq_len(nl)); hc$labels <- ctn
tree <- hclust_to_tree(hc, ctn)

n <- 1500L
leaf_ll <- matrix(rnorm(n * tree$n_leaf, -3, 0.6), n, tree$n_leaf)
ll  <- vapply(tree$leaf_all, function(v) rowSums(leaf_ll[, v, drop = FALSE]), numeric(n)) +
       matrix(rnorm(n * tree$n_col, 0, 0.3), n, tree$n_col)
llt <- ll - abs(matrix(rnorm(n * tree$n_col, 0, 0.15), n, tree$n_col))
fit <- estimate_tree_prior(tree, ll, llt, tol = 1e-8, max_iter = 120)

mu_full <- matrix(rnorm(n * tree$n_col, 0, 1.8), n, tree$n_col)
sigma_full <- matrix(runif(n * tree$n_col, 0.05, 0.25), n, tree$n_col)
Nmean <- rep(40, n)

m <- yaml::read_yaml(cft_default_config())$markers
m$min_read_depth <- 5
m$lowest_p <- 0.5
run <- function(mm) select_markers(tree, mm, ll, llt, mu_full, sigma_full,
                                   fit$pi.j, fit$pi.t, Nmean, seq_len(n))
key <- function(d) sort(paste(d$index, d$node, sep = "|"))

d <- run(m)
semi <- d[d$onevsrest == 0L, , drop = FALSE]
stopifnot(all(is.na(d$semi_pair[d$onevsrest == 1L])), !anyNA(semi$semi_pair))
cat(sprintf("\n  semi rows %d: pair %d, single %d\n",
            nrow(semi), sum(semi$semi_pair), sum(!semi$semi_pair)))
stopifnot(any(semi$semi_pair), any(!semi$semi_pair))

pa  <- calculate.p.gj(tree$children, fit$pi.j, ll, tree$layer)
p.t <- calculate.p.t(fit$pi.t, fit$pi.j, tree, ll, llt, pa[[2]])
both_blocks <- function(g, v) {
  par <- tree$parent[v]; sib <- setdiff(tree$children[[par]], v)[1]
  pr <- pa[[1]][g, v] * pa[[1]][g, sib] * p.t[g, 1]
  a <- par
  repeat { pr <- pr * (1 - pa[[1]][g, a]); if (a == tree$n_node) break; a <- tree$parent[a] }
  pr
}
ref <- mapply(both_blocks, semi$local_idx, semi$node) >= m$lowest_p
stopifnot(identical(unname(ref), semi$semi_pair))
cat("  label matches P(node and sibling both blocks) >= lowest_p   PASS\n")

m2 <- m; m2$gate_requires_sibling_block <- TRUE
d2 <- run(m2)
semi2 <- d2[d2$onevsrest == 0L, , drop = FALSE]
stopifnot(all(semi2$semi_pair), identical(key(semi2), key(semi[semi$semi_pair, ])))
cat("  semi_pair rows == rows selected with gate_requires_sibling_block   PASS\n")
cat("\nSEMI-PAIR TESTS PASSED\n")
