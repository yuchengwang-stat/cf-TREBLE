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

m <- yaml::read_yaml(cft_default_config())$markers
m$min_read_depth <- 5
m$lowest_p <- 0.5
run <- function(mm) select_markers(tree, mm, ll, llt, mu_full, sigma_full,
                                   fit$pi.j, fit$pi.t, rep(40, n), seq_len(n))

stopifnot(!length(m$strict_celltypes))
base <- run(m)
m_hi <- m; m_hi$strict_rel_distance <- 99
stopifnot(identical(run(m_hi), base))
cat("\n  empty strict_celltypes: strict_rel_distance has no effect   PASS\n")

m_idx <- m; m_idx$strict_celltypes <- list(3L, 7L)
m_nm  <- m; m_nm$strict_celltypes  <- list(ctn[3], ctn[7])
a <- run(m_idx); b <- run(m_nm)
stopifnot(identical(a, b), !identical(a, base))
cat(sprintf("  names and indices select the same markers (%d rows)       PASS\n", nrow(a)))

m_bad <- m; m_bad$strict_celltypes <- list("NotACellType")
stopifnot(inherits(try(run(m_bad), silent = TRUE), "try-error"))
cat("  an unknown cell type is rejected                           PASS\n")
cat("\nSTRICT-CELLTYPES TESTS PASSED\n")
