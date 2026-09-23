# cfTREBLE

Tree-based Bayesian estimation of cell-type methylation signatures, and
deconvolution of cell-free DNA against them.

A reference panel of bisulfite-sequenced sorted cell types goes in; a set of
marker CpGs and a deconvolution signature come out, and any `.beta` sample can
then be resolved into cell-type fractions. The work is split into seven stages
driven by one configuration file, each restartable, designed to run as a SLURM
job array.

## Installation

R >= 4.1, and:

```r
install.packages(c("data.table", "pbv", "pracma", "RhpcBLASctl",
                   "mixsqp", "yaml", "jsonlite"))
# these are all of them -- nothing else is loaded, installed or from a checkout
BiocManager::install("impute")   # optional, only for `tree --from bed`
```

Then:

```r
remotes::install_github("yuchengwang-stat/cfTREBLE")
# or, from a checkout:
R CMD INSTALL cfTREBLE
```

The command-line driver ships inside the package:

```bash
CFT=$(Rscript -e 'cat(system.file("scripts", "cftreble", package = "cfTREBLE"))')
"$CFT" --help
```

Put it on your `PATH` as `cftreble` if you like; the rest of this file assumes
that. The scripts also run from a git checkout without installing.

## Input data

**Genome build: hg38.** A CpG is identified throughout by its row number in a
genome-wide CpG list, ordered by chromosome and coordinate, and by nothing else.
That list is the hg38 one used by
[wgbs_tools](https://github.com/nloyfer/wgbs_tools). The reference panel, the
signature and every sample being deconvolved must be indexed against the same
list; a panel on hg19, or indexed some other way, has to be re-indexed first.
Nothing in the pipeline can detect a mismatch, except an index that runs off the
end of a `.beta` file.

**Autosomes only.** The list used here covers chr1 to chr22 and excludes chrX,
chrY and chrM, which comes to 27,852,739 CpGs -- the number the shipped
`data.chunks` and `chunk_stride` describe. This matters more than it looks:
a list that includes the sex chromosomes is a different list, and every index
after the first extra CpG shifts, silently. Whatever you choose, build the
reference panel and the `.beta` files you deconvolve the same way, and check the
sizes line up -- a `.beta` is exactly 2 bytes per CpG, so an autosome-only hg38
file is 55,705,478 bytes.

Four things go in, all pointed at from the configuration file.

**Reference counts.** `M` (methylated) and `N` (total) integer matrices, rows =
CpG, columns = reference sample, saved as `.rds` and split into chunks of equal
height. A CpG's row number within the whole genome-wide list is its identifier
and must match the row number in the `.beta` files you later deconvolve.

**Sample list.** A named list, one entry per cell type, holding that cell type's
sample names. The column order of `M` and `N` must equal `unlist(sample_list)`.

**Priors.** `mu`, a 3-component normal mixture with `$mu`, `$sigma`, `$lambda`;
and `sigma`, a grid with `$grid` and `$pi`.

**Cell-type tree.** Either supplied as `children` and `parent` lists, or built
from the reference counts with `cftreble tree --from ref`.

**Samples to deconvolve** are `.beta` files, the format
[wgbs_tools](https://github.com/nloyfer/wgbs_tools) writes: a flat stream of
`uint8` (methylated, total) pairs, one pair per CpG, one row per CpG of the
genome-wide list, with no header
([format description](https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md)).
The CpG count is read from the file size, so nothing has to be configured to
match.

To see every one of these objects as a concrete file, run
`cftreble example --out demo`. It writes a small but complete dataset -- sample
list, tree, priors, count matrices, a `.beta` sample with its true composition --
plus a configuration that runs on them. Copy the shapes from there.

## Quick start

```bash
cftreble demo
```

builds that dataset, runs all seven stages on it, and prints the estimated
cell-type fractions next to the proportions the sample was mixed at. It echoes
each command as it runs, so the transcript is the sequence to copy. Takes about
a minute and needs no data of your own.

On your own data:

```bash
cp $(Rscript -e 'cat(system.file("config", "default.yaml", package = "cfTREBLE"))') my.yaml
# edit the paths at the top of my.yaml

cftreble --config my.yaml check          # validate the inputs first
cftreble --config my.yaml submit         # stages 1-6 as a SLURM chain
cftreble --config my.yaml status         # what has finished
```

`submit` chains stages 1 to 6 with SLURM dependencies and turns the per-chunk
ones into job arrays, using the resources in the `slurm:` block of the config.
Stage 7 is not chained -- it needs the samples you want deconvolved, so run it
yourself once the signature exists. To run any stage by hand:

```bash
cftreble --config my.yaml tree
cftreble --config my.yaml likelihood --chunk 3
cftreble --config my.yaml deconvolve --beta '/path/to/*.beta' --out fractions.csv
```

## Stages

| # | stage | input -> output | parallel |
|---|-------|-----------------|----------|
| 1 | `tree` | builds and validates the cell-type tree every later stage reads | - |
| 2 | `likelihood` | per-node marginal log-likelihood for one CpG chunk | array over chunks |
| 3 | `prior` | EM for the Bernoulli and categorical priors on the tree | - |
| 4 | `markers` | marker CpG selection for one chunk | array over chunks |
| 5 | `signature` | merges chunks, applies the per-label caps, writes the signature | - |
| 6 | `subject` | per-reference-subject signature over the marker CpGs | array over chunks |
| 7 | `deconvolve` | `.beta` sample(s) -> cell-type fractions | - |

Every stage writes a `manifest.json` and skips work it has already done, so a
job that hits the wall clock can be resubmitted as-is.

## What it costs

Measured on Bridges-2 `RM-shared`, 26 cores, a 207-sample reference over a
187-column tree. The likelihood stage is essentially the whole bill; everything
after it is minutes.

| | |
|---|---|
| likelihood | 560 CpG/min on 26 cores |
| one chunk of 900,000 CpGs | ~27 hours |
| all 31 chunks as a job array | ~27 hours wall, **~21,600 core-hours** |
| prior, markers, signature | minutes |
| subject, per chunk | ~45 min |
| deconvolve, per sample | seconds |

Three things follow.

- 21,600 core-hours is a real bite out of an allocation. Decide before
  submitting, and consider whether the marker scan needs all 27.8M CpGs or
  whether a read-depth or variability pre-filter can cut the input first.
- `slurm.time` defaults to `48:00:00` because a chunk needs ~27 h. The stage is
  resumable per batch, so a job that hits the wall can simply be resubmitted --
  finished batches are skipped.
- Memory peaked at 8.6 GB. The dominant term is the chunk itself (900,000 x 207
  integers, ~745 MB each for `M` and `N`), not the batch size, so a full-size
  run is not much larger. `--ntasks-per-node=26` buys ~52 GB on `RM-shared`.

Scale these by your own panel: the cost goes with CpGs x tree columns x samples.

## The tree

One object, written by stage 1, that everything downstream reads:

```
celltypes  children  parent  leaf  leaf_all  comp_of  layer
n_leaf     n_node    n_col
```

`layer`, the complement columns and `n_col` are derived from `children`, not
configured. Complements carrying no information are dropped: the root's
complement is empty, and the root's two children are each other's complement. A
48-leaf tree gives 95 nodes and 187 likelihood columns. Swapping in a different
tree needs no change downstream.

## Marker types

Every marker is exactly one of three kinds, reported in the `kind` column.

| kind | meaning |
|---|---|
| `celltype` | one cell type against all the others |
| `class` | one internal node -- a group of cell types -- against all the others |
| `semi` | neither: the CpG splits the tree into blocks, without any single block standing against everything else |

The `target` column names what the marker is a marker *for*: a cell type by
name, an internal node as `node<k>`.

The same CpG can be a marker for more than one target, and `target` then lists
them separated by `;`, with `n_target` counting them.

```
index   target                  kind      n_target
1234    CellType03              celltype  1
5678    node70                  class     1
8116    CellType34;CellType35   semi      2
```

`celltype` and `class` markers always have `n_target` 1, since the one-vs-rest
test picks a single winning partition per CpG.

Stage 5 writes `signature.txt.gz` with all of them, plus
`signature_one_vs_rest.txt.gz` (`celltype` and `class`) and
`signature_semi.txt.gz` separately. Each is a table of the metadata columns
above followed by one beta value per cell type.

## Deconvolution

Stage 7 fits BLEND to each `.beta` sample over the signature CpGs, with
`cf_BLEND()`. A cell-free sample is modelled as a mixture over cell types, each
of which is itself a mixture over the reference subjects of that type. The
methylated and unmethylated references are stacked into one convex problem and
solved with sequential quadratic programming, giving a cellular fraction per cell
type and a mixing proportion over that type's reference panel. Stage 7 writes the
cellular fractions; `cf_BLEND()` returns both.

On data simulated from the model it recovers the mixing proportions to about
2e-03 (`inst/tests/test_blend.R`).

### Per-subject references

BLEND takes `phi` as one matrix per cell type, whose columns are that cell type's
reference subjects -- if `Small-Int-Ep` has four samples, its entry is a G x 4
matrix. By default stage 7 gives each cell type a single column, the signature's
own value for it, which is a reference panel of one subject per type.

Stage 6 computes the real thing: a posterior methylation estimate per reference
subject at each marker CpG. Pass `--subject` to deconvolve against it.

```bash
cftreble --config my.yaml deconvolve --beta sample.beta --subject
cftreble --config my.yaml deconvolve --beta sample.beta --subject /other/run/subject
```

With it, `cf_BLEND()`'s `ref_mixing_prop` becomes a real quantity -- how the
sample's signal distributes over that cell type's subjects. Without it the
columns within a cell type are identical, so the mixing proportions are an
arbitrary split of a tie, while the cellular fractions are unaffected.

A reference subject with no reads at a marker CpG has no estimate there.
`--na drop` (the default) drops those CpGs; `--na celltype` substitutes the cell
type's own value for that subject at that CpG.

Whether this helps is an empirical question and depends on how much the subjects
of a cell type actually differ. On the synthetic demo, where they differ only by
noise around a common truth, it is slightly worse (RMSE 0.018 against 0.016).

### Using a signature built elsewhere

No reference data ships with this package, so there is no signature to hand out.
Deconvolution does not need one from this pipeline: `--signature` also accepts a
plain table, and nothing else about the reference panel is required.

The table needs an `index` column giving each CpG's row number in the `.beta`
files you are deconvolving, and one column per cell type holding that cell
type's beta value in `[0, 1]`. Tsv or csv, optionally gzipped. Any other column
this pipeline writes (`target`, `kind`, `n_target`, `mu_*`, `sigma_*`) is
ignored, so `signature.txt.gz` can be read straight back in.

```
index   CellType01  CellType02  CellType03
7       0.1137      0.1137      0.8879
2       0.3830      0.6880      0.6895
```

```bash
cftreble --config my.yaml deconvolve \
  --signature my_signature.tsv \
  --beta '/path/to/*.beta' \
  --out fractions.csv

# restrict the reference to some of its cell types
cftreble --config my.yaml deconvolve \
  --signature my_signature.tsv \
  --celltypes Monocyte,Neutrophil,Hepatocyte \
  --beta sample.beta
```

The only thing that has to line up is the CpG indexing: `index` is a row number
into the `.beta` file, so the signature and the samples must be on the same CpG
list. An index past the end of a `.beta` is an error rather than a silent wrap.

One thing not to worry about: the signature is per cell type, and giving a cell
type several reference columns instead of one leaves the cellular fractions
unchanged -- it only affects the mixing proportions within a type
(`inst/tests/test_deconvolve.R` checks this to 1e-10).

## Configuration

One file; copy `inst/config/default.yaml` and pass it with `--config`. `${key}`
expands against the top-level scalars, so paths are written once.

Notes on individual settings:

- `markers.max_per_celltype` and `max_per_node` cap how many markers one label
  may contribute, applied after the chunks are merged. Marker counts are very
  uneven across cell types, and without a cap two of them can hold most of the
  one-vs-rest markers.
- `markers.sibling_contrast` is an additional, off-by-default scan that compares
  the two children of each internal node directly; `background` sets what must
  hold of the remaining cell types.
- `priors.sigma_truncate: auto` selects the grid points below 0.5.
- `likelihood.ncores` should match `slurm.ntasks`. `block_size` is capped at
  `ceil(rows / ncores)` so a short batch still uses every core.
- `subject.weight_coverage` (default 0.9999) prunes the partition sum. Every
  value being averaged is a probability, so the error is bounded by
  `1 - coverage`; measured worst case 2.9e-05. Set `1.0` to compute it exactly.
- `prior_em` runs over every CpG, accumulating chunk by chunk rather than
  holding them all in memory. `pi.t` has one component per partition and the
  rare ones are what marker selection keys on, so subsampling would be noisy
  where it matters; the streaming form is identical to the batch computation,
  not an approximation (`inst/tests/test_em_streaming.R`).

## Testing

```bash
cftreble demo                              # all seven stages on a generated dataset
cftreble selftest                          # synthetic, end to end
cftreble validate --cpgs 3000              # synthetic with known truth
Rscript inst/tests/test_em_streaming.R     # streaming EM == batch EM
Rscript inst/tests/test_subject.R          # subject weights are a convex combination
Rscript inst/tests/test_blend.R            # cf_BLEND recovers a known mixture
Rscript inst/tests/test_deconvolve.R       # deconvolution against an external signature
Rscript inst/tests/test_submit.R           # submit builds the right SLURM chain
Rscript inst/tests/test_sibling_contrast.R # sibling_contrast: inert off, additive on
```

`validate` plants methylation truth on the tree so all three kinds of marker
exist, runs the pipeline, deconvolves mixtures at known proportions, and reports
marker precision and recall against what it planted alongside the deconvolution
error. It exercises the implementation under the model's own assumptions and
says nothing about real WGBS artefacts.

## Data availability

The reference panel behind the published signature combines sources that are not
all redistributable, so no reference data is included here and there is no script
that regenerates that signature. `cftreble demo` and `cftreble validate` run the
whole pipeline on data they generate themselves.

Every threshold used in selection is a setting in the configuration file rather
than a constant in the code. The shipped values are defaults; a different panel,
sequencing depth or tree will want different ones.

## Earlier implementations

Two directories under `inst/` keep the code the method was developed with. The
package neither loads nor compiles either of them; each has a README of its own.

`inst/legacy-R/` is the `R/` directory of the first implementation, verbatim.
Four of its ten files were carried into the package unchanged and only renamed;
the other six were reworked into `R/tree.R`, `R/likelihood.R`, `R/prior_em.R`,
`R/markers.R` and `R/deconvolve.R`. Its README gives the file-by-file mapping.

`inst/legacy-cpp/` holds the C++: `fragment.cpp`, which computes read-level
likelihoods without collapsing reads to per-CpG counts, and
`our_celfie_rcpp.cpp`, an Armadillo rewrite of the CelFiE EM. They are not in
`src/` because they need `RcppArmadillo`, `pbv` and `BH`, which the package does
not depend on; their README shows how to compile them on their own.

## License

Not yet chosen -- see `LICENSE`. Patent matters relating to the method are still
open, so the default of copyright applies for now.
