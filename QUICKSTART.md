# Quickstart

Every command below was run on a cluster from a fresh clone of this repository;
the output shown is what it actually printed.

## 0. If you have never seen these inputs before

Do not guess the formats from prose. Write a complete, valid, miniature example
and look at it:

```bash
CFT=$(Rscript -e 'cat(system.file("scripts", "cftreble", package = "cfTREBLE"))')
"$CFT" example --out demo
```

```
Wrote a complete example under demo

  reference/sample_list.rds        named list, 8 cell types x 3 samples
  reference/tree_children.rds      15 nodes; leaves 1..8 are c(i), internal are c(left, right)
  reference/tree_parent.rds        15 integers; parent[15] = 15 (the root is its own parent)
  reference/mu_prior.rds           list($mu, $sigma, $lambda), 3 components
  reference/sigma_prior.rds        list($grid, $pi), grid ascending and starting at 0
  reference/chunks/{M,N}_k.rds     2 chunks of 400 CpGs x 24 samples (integer counts)
  reference/plasma_example.beta    uint8 (methylated, total) pairs, one per CpG
  reference/plasma_example_truth.rds  the proportions it was mixed at
  example.yaml                     a config pointing at all of it
```

Every object the pipeline needs, small enough to open:

```bash
Rscript -e 'str(readRDS("demo/reference/tree_children.rds"))'
```

The whole pipeline runs on it in a couple of minutes, and the plasma sample's
true composition is stored next to it, so you can see the answer come back:

```
  cell type         truth      est
  CellType01       0.5000   0.4858
  CellType02       0.2000   0.2044

  RMSE 0.0161   Pearson 0.994
```

## 1. Install

```bash
git clone https://github.com/yuchengwang-stat/cf-TREBLE.git
R CMD INSTALL --no-docs -l /path/to/your/Rlib cf-TREBLE
```

or, in R:

```r
remotes::install_github("yuchengwang-stat/cf-TREBLE")
```

On a cluster, point `R_LIBS_USER` at the library you installed into, in every
batch script:

```bash
export R_LIBS_USER="/path/to/your/Rlib:$R_LIBS_USER"
```

## 2. Find the driver

```bash
CFT=$(Rscript -e 'cat(system.file("scripts", "cftreble", package = "cfTREBLE"))')
"$CFT"
```

```
cf-TREBLE pipeline driver.

  cftreble <stage> [--config FILE] [stage args]      run a stage here and now
  cftreble submit  [--config FILE] [--from STAGE]    submit the whole chain to SLURM
  cftreble smoke   [--rows N] [--beta FILE]          tiny end-to-end run, no SLURM
  cftreble status  [--config FILE]                   what has finished so far
  cftreble selftest [--cpgs N] [--tree-lists DIR]    synthetic end-to-end check
  cftreble validate [--cpgs N] [--plasma-depth D]    accuracy vs known truth

Stages: tree | likelihood | prior | markers | signature | subject | deconvolve
```

## 3. Check the install before touching your data

```bash
"$CFT" selftest --cpgs 60
```

```
======== checks ========
  tree validates                               PASS
  n_col == n_node + #complements               PASS
  layer partitions internal nodes              PASS
  EM produced finite pi                        PASS
  pi.t sums to 1                               PASS (1.000000)
  signature non-empty                          PASS (30 CpGs)
  signature beta in [0,1]                      PASS
  signature has one column per cell type       PASS
  fractions sum to 1                           PASS (1.000000)
  fractions non-negative                       PASS
  dominant cell type recovered                 PASS (got CT1, truth CT1)

 SELFTEST PASSED
```

Ten seconds, synthetic data, no cluster. It builds a small tree, runs every
stage, and deconvolves a mixture whose composition it knows. If this fails, the
environment is wrong and nothing else is worth trying.

## 4. Point a config at your data

```bash
cp $(Rscript -e 'cat(system.file("config", "default.yaml", package = "cfTREBLE"))') my.yaml
```

Two lines decide everything else:

```yaml
base:      /your/output/dir        # where results go
readonly:  /your/reference/dir     # inputs, never written to
```

The rest of the paths are written as `${readonly}/...`, so those two roots are
the only things most sites change.  `cftreble example` writes a complete set of
inputs and a matching config if you would rather start from a working one.

What has to exist under `readonly`:

| key | what it is |
|---|---|
| `data.m_chunks` / `data.n_chunks` | `sprintf` patterns for the methylated and total count matrices, one file per chunk, rows = CpG, columns = reference sample |
| `data.sample_list` | named list, one entry per cell type, holding its sample names -- **column order of M/N must equal `unlist(sample_list)`** |
| `priors.mu` | 3-component normal mixture with `$mu`, `$sigma`, `$lambda` |
| `priors.sigma` | grid with `$grid` and `$pi` |
| `tree.from_lists` | `children` and `parent` lists -- or skip them and use `tree --from ref` |

## 4b. Check the inputs before spending anything

```bash
"$CFT" --config my.yaml check
```

```
  ok   sample_list is a named list                8 entries
  ok   leaves are nodes 1..n_leaf
  ok   tree is binary                             15 nodes, 8 leaves; expected 15 nodes
  ok   one leaf per cell type                     8 leaves, 8 cell types
  ok   root is its own parent                     parent[15] = 15
  ok   tree builds and validates                  8 leaves, 15 nodes, 27 likelihood columns, 5 layers
  ok   mu prior weights sum to 1                  sum = 1.000000
  ok   grid is ascending                          0 0.1 0.15 0.2 0.25 0.35 0.5 0.6 1
  ok   chunk 1: columns match sample_list         24 columns, sample_list implies 24
  ok   chunk 1: M <= N everywhere
  ...
  26 checks, 0 failed
```

Twenty-six checks, seconds, before the first expensive job. It catches the
mistakes that otherwise surface late or not at all -- a tree whose leaves are not
the first `n_leaf` nodes still "works" and quietly produces a wrong signature.
`--deep` opens every chunk instead of the first.

## 5. Build the tree

```bash
"$CFT" --config my.yaml tree
```

```
tree: 48 leaves, 95 nodes, 187 likelihood columns, 8 layers
```

Seconds. `layer`, the informative complement set and `n_col` are all derived
from `children` -- nothing to configure. If the tree does not validate, this
stops here rather than producing something subtly wrong later.

No tree to hand it? Build one from the reference counts:

```bash
"$CFT" --config my.yaml tree --from ref
```

## 6. Run the rest

Whole chain as a SLURM dependency chain:

```bash
"$CFT" --config my.yaml submit
"$CFT" --config my.yaml status
```

```
  tree         done
  likelihood   done (2/2 chunks)
  prior        done
  markers      done (2/2 chunks)
  signature    done
  subject      done (2/2 chunks)
  deconvolve   done
```

Array stages report chunks finished rather than a single flag, so a run that is
half-done looks half-done.

One stage at a time, if you would rather watch:

```bash
"$CFT" --config my.yaml likelihood --chunk 3
"$CFT" --config my.yaml prior
"$CFT" --config my.yaml markers --chunk 3
"$CFT" --config my.yaml signature
"$CFT" --config my.yaml subject --chunk 3
"$CFT" --config my.yaml deconvolve --beta '/path/to/*.beta' --out fractions.csv
```

## 7. Try it small first

```bash
"$CFT" --config my.yaml smoke --rows 3000
```

Runs every stage on 3,000 real CpGs without SLURM. Use `--sample-rows` rather
than `--max-rows` for anything where the answer matters: `--max-rows` takes the
head of a chunk, which is one contiguous stretch of one chromosome, while markers
are scattered genome-wide.

## What it costs

The likelihood stage is essentially the whole budget; everything after it takes
minutes. The chunks run as an array, so wall time is one chunk rather than all
of them.

## Resuming

Stages are resumable at batch granularity. A rerun checks what is already on
disk -- including whether it covers the rows this run wants -- and skips only what
genuinely matches. Killing a job and resubmitting the same command is safe.
