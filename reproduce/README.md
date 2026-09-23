# Reproducing the deconvolution results

Scripts used for the paper, kept apart from the package: neither
`R CMD INSTALL` nor `remotes::install_github` installs this folder.

## Deconvolution demo

`deconvolution_demo.R` deconvolves an example cfDNA sample with the published
subject-specific signature (628,220 marker CpGs) using `cf_BLEND`. The signature
covers 46 of the paper's 47 cell types and their 196 reference samples; the
chorionic villus (CVS) samples are not publicly released, so CVS is omitted.
From the repository root:

```bash
Rscript reproduce/deconvolution_demo.R
```

On first use it downloads two files attached to the
[v1.0.0 release](https://github.com/yuchengwang-stat/cf-TREBLE/releases/tag/v1.0.0)
into `reproduce/data/`:

| file | size | contents |
|---|---|---|
| `cfTREBLE_signature_630k_subjspec.txt.gz` | 250 MB | CpG index and hg38 coordinates, then one column per reference sample |
| `demo_data.beta` | 59 MB | example cfDNA sample in [wgbs_tools](https://github.com/nloyfer/wgbs_tools) `.beta` format |

`data/sample_list.RDS` maps each cell type to its reference samples, in the
column order of the signature.

The estimated fractions are written to `reproduce/output/demo_result.csv` and
compared with `reproduce/demo_result.csv`. The run takes about 11 minutes and
7 GB of memory on a laptop.

## Simulated cfDNA mixtures

The simulations mix fragments from held-out reference samples (Loyfer et al.,
2023). Two scripts generate them.

1. `generate_pseudo_plasma.R` draws one set of ground-truth proportions --
   fixed proportions for T, NK and B cells, monocytes, macrophages and
   granulocytes, plus eight other cell types chosen at random -- and samples
   `nfrag` fragments per mixture from each cell type's `.pat.gz` file, writing
   20 replicates. The proportions are saved as `true_prop_<param>.RDS`.

   ```bash
   Rscript reproduce/generate_pseudo_plasma.R test_samples.rds pat_dir out_dir <param> <nfrag> [seed]
   ```

   `test_samples.rds` is a list of length 48, one entry per cell type in tree
   order, holding the name of the held-out sample or `NULL`. It lists the 48
   reference cell types before monocytes and macrophages are combined into the
   single cell type used in the paper.

2. `generate_beta.sh` merges the 14 fragment files of one replicate into a
   sorted `.pat` file, keeps a random half of its reads, and converts it to
   `.beta` with wgbs_tools.

   ```bash
   bash reproduce/generate_beta.sh out_dir <param> <replicate>
   ```

The resulting `.beta` files are deconvolved as in `deconvolution_demo.R`.
