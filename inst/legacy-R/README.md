# Original R sources

The `R/` directory of the first cfTREBLE implementation, kept verbatim as a
record of the work. Nothing here is loaded: these files live under `inst/`, not
`R/`, so the package never sources them and the names below do not collide with
the installed ones.

The current package covers the same ground with a different file layout. Four
files were carried across unchanged, only renamed; the other six were reworked.

## Carried across verbatim

| Original | Now |
| --- | --- |
| `FastDerivatives_Post.R` | `R/math_fastderivatives_post.R` |
| `FastDerivatives_Post_UpperNodes.R` | `R/math_fastderivatives_post_uppernodes.R` |
| `PosteriorExpectations.R` | `R/math_posteriorexpectations.R` |
| `PosteriorExpectations_IndSig.R` | `R/math_posteriorexpectations_indsig.R` |

These four are byte-identical to the installed versions. They are duplicated
here so the original directory is complete rather than partial.

## Reworked

| Original | What it held | Now |
| --- | --- | --- |
| `BuildClusterTree.R` | `BuildClusterTree` -- hierarchical clustering of the reference panel into the cell-type tree | `R/tree.R` |
| `CalculateLikelihood.R` | `CalculateLikelihood`, `CalculateIndp` -- per-node marginal likelihoods over the panel | `R/likelihood.R` |
| `EstimateTreePrior.R` | `append_complements`, `EstimateTreePrior` -- Bernoulli and categorical priors on the tree | `R/prior_em.R` |
| `TreeBasedEM.R` | `calculate.p.gj`, `calculate.w.gj`, `calculate.p.t`, `calculate.pi.j` -- the EM updates behind that estimation | `R/prior_em.R` |
| `SelectMarkerCpGs.R` | `SelectMarkerCpGs` -- marker CpG selection from the fitted posteriors | `R/markers.R` |
| `cf_BLEND.R` | `cf_BLEND`, `cf_BLEND_ONE` -- the deconvolution itself | `R/deconvolve.R` |

`cf_BLEND` is the one reworked file whose body did not change: the installed
copy differs only in replacing `rlist::list.cbind(.)` with
`do.call(cbind, .)`, which removed a dependency without changing the result.
