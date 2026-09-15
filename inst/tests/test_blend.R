#!/usr/bin/env Rscript
## Both BLEND solvers fit the same mixture model, so on data simulated from that
## model they must recover the same cellular fractions and the same truth.
suppressMessages(library(cfTREBLE))
set.seed(7)

G  <- 4000L                      # signature CpGs
M  <- c(A = 5L, B = 4L, C = 6L)  # reference subjects per cell type
phi <- lapply(M, function(m) matrix(rbeta(G * m, 0.6, 0.6), G, m))

truth  <- c(A = 0.55, B = 0.30, C = 0.15)
mixing <- lapply(M, function(m) { w <- runif(m); w / sum(w) })
p  <- Reduce(`+`, Map(function(t, f, w) f * (phi[[t]] %*% w), names(M), truth, mixing))
D  <- rep(400L, G)
X  <- rbinom(G, D, pmin(pmax(p, 1e-6), 1 - 1e-6))

fit <- cf_BLEND(X, D, phi)

cat(sprintf("  %-4s truth %.4f   estimated %.4f\n",
            names(truth), truth, fit$cellular_frac), sep = "")
cat(sprintf("  max |est - truth| = %.2e ; fractions sum to %.10f\n",
            max(abs(fit$cellular_frac - truth)), sum(fit$cellular_frac)))

stopifnot(max(abs(fit$cellular_frac - truth)) < 0.02)
stopifnot(abs(sum(fit$cellular_frac) - 1) < 1e-8)
stopifnot(identical(names(fit$cellular_frac), names(phi)))
## the mixing proportions are per cell type, one weight per reference column
stopifnot(vapply(fit$ref_mixing_prop, length, 0L) == M,
          all(abs(vapply(fit$ref_mixing_prop, sum, 0) - 1) < 1e-8))

## Equal-sized reference panels take a different path through the internal
## sapply(), which simplifies to a matrix when every cell type has the same
## number of reference columns.
Me  <- c(A = 3L, B = 3L, C = 3L)
phe <- lapply(Me, function(m) matrix(rbeta(G * m, 0.6, 0.6), G, m))
mxe <- lapply(Me, function(m) { w <- runif(m); w / sum(w) })
pe  <- Reduce(`+`, Map(function(t, f, w) f * (phe[[t]] %*% w), names(Me), truth, mxe))
Xe  <- rbinom(G, D, pmin(pmax(pe, 1e-6), 1 - 1e-6))
eqe <- cf_BLEND(Xe, D, phe)
cat(sprintf("  equal panels: max |est - truth| = %.2e\n",
            max(abs(eqe$cellular_frac - truth))))
stopifnot(vapply(eqe$ref_mixing_prop, length, 0L) == Me,
          max(abs(eqe$cellular_frac - truth)) < 0.02)
cat("BLEND TESTS PASSED\n")
