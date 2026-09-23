#!/usr/bin/env Rscript
## Summarise what `cftreble demo` produced, against the truth it was built from.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a   <- args_parse(list(out = ""))
if (!nzchar(a$out)) die("--out <demo directory> is required")
cfg <- load_config(file.path(a$out, "example.yaml"))
S   <- readRDS(file.path(run_dir(cfg, "signature", FALSE), "signature.rds"))

cat(sprintf("  signature   %d CpGs x %d cell types\n", nrow(S$table), ncol(S$beta)))
for (k in c("celltype", "class", "semi_pair", "semi_single"))
  cat(sprintf("    %-11s %5d CpGs\n", k, sum(grepl(k, S$table$kind, fixed = TRUE))))

est <- unlist(read.csv(file.path(run_dir(cfg, "deconvolve", FALSE),
                                 "cell_fractions.csv"), row.names = 1))
tru <- readRDS(file.path(a$out, "reference", "plasma_example_truth.rds"))[names(est)]
cat("\n  cell type        mixed at   estimated\n")
cat(sprintf("  %-14s %8.3f %11.3f\n", names(tru), tru, est), sep = "")
cat(sprintf("\n  RMSE %.4f   max abs error %.4f   Pearson %.4f\n",
            sqrt(mean((tru - est)^2)), max(abs(tru - est)), cor(tru, est)))
