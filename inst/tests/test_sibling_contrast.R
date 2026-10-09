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
run <- function(mm) select_markers(tree, mm, ll, llt, mu_full, sigma_full,
                                   fit$pi.j, fit$pi.t, Nmean, seq_len(n))
key <- function(d) if (is.null(d)) character(0) else
  sort(paste(d$index, d$node, d$onevsrest, sep = "|"))

off <- run(m); ko <- key(off)
stopifnot(!isTRUE(m$sibling_contrast$enabled))

m_none <- m; m_none$sibling_contrast <- NULL
stopifnot(identical(run(m_none), off))
cat(sprintf("\n  contrast off == block absent           PASS  (%d rows)\n", length(ko)))

any_added <- FALSE
for (bg in c("iqr", "joint_p", "span", "both", "none")) {
  mm <- m; mm$sibling_contrast$enabled <- TRUE; mm$sibling_contrast$background <- bg
  kn <- key(run(mm))
  added <- length(setdiff(kn, ko)); lost <- length(setdiff(ko, kn))
  cat(sprintf("  background %-8s %4d -> %4d rows   added %4d  lost %d   %s\n",
              bg, length(ko), length(kn), added, lost,
              if (lost == 0L && !anyDuplicated(kn)) "PASS" else "FAIL"))
  stopifnot(lost == 0L, !anyDuplicated(kn))
  any_added <- any_added || added > 0L
}
stopifnot(any_added)

mm <- m; mm$sibling_contrast$enabled <- TRUE; mm$sibling_contrast$background <- "nonsense"
stopifnot(inherits(try(run(mm), silent = TRUE), "try-error"))
cat("  an unknown background is rejected     PASS\n")
cat("\nSIBLING-CONTRAST TESTS PASSED\n")
