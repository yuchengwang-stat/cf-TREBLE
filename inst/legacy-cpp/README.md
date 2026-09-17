# Legacy C++ sources

Earlier C++ implementations written during the development of cfTREBLE. They are
kept here as a record of the work; **no function in the installed package calls
them**, and they are deliberately not in `src/`, so `R CMD INSTALL` does not try
to compile them.

They are not compiled with the package because they need headers the package
itself does not depend on (`RcppArmadillo`, `pbv`, `BH`). To use one, compile it
on its own after installing those:

```r
install.packages(c("RcppArmadillo", "pbv", "BH"))
Rcpp::sourceCpp(system.file("legacy-cpp", "fragment.cpp", package = "cfTREBLE"))
```

## `fragment.cpp`

Fragment-level likelihood machinery: reads are kept whole rather than collapsed
to per-CpG counts, so the correlation between CpGs carried on the same read is
retained. 1,652 lines, ten exported entry points.

| Function | What it does |
| --- | --- |
| `log_lik_laplace_cpp` | Laplace-approximated log-likelihood over a tree of cell types |
| `estimate_para` | Parameter estimation on that tree |
| `find_intervals_cpp` | Split a chromosome into intervals by CpG density |
| `find_pairs` / `find_pairs_subj` / `find_pairs_cpp` | Co-occurring CpG pairs within reads |
| `find_overlap_fragments` | Reads overlapping a given interval |
| `findIndicesInRange` | Binary search over a sorted position vector |
| `calculate_1_fragment_lik_cpp` | Likelihood of one read under a correlation model |
| `calculate_fragment_lik_cpp` | The same across all reads of an interval |

Dependencies: `RcppArmadillo`, `pbv` (bivariate normal CDF), `BH` (Boost, for
the incomplete beta function).

## `our_celfie_rcpp.cpp`

An Armadillo rewrite of the CelFiE EM algorithm: alternates between cell-type
proportions and reference methylation rates until the proportions stop moving.
Exports `my_celfie_rcpp(X, D_X, W, D_W, n_iter, thres)`.

Dependencies: `RcppArmadillo`.
