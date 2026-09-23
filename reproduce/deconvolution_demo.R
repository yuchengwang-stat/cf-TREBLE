#!/usr/bin/env Rscript
## Deconvolve the example cfDNA sample with the published subject-specific
## signature (about 630k marker CpGs). Run from the repository root:
##
##   Rscript reproduce/deconvolution_demo.R
##
## The signature and the example .beta file are downloaded from the v1.0.0
## release on first use. Estimated fractions are written to
## reproduce/output/demo_result.csv and compared with reproduce/demo_result.csv.
suppressMessages({ library(data.table); library(cfTREBLE) })

dir <- "reproduce/data"
out <- "reproduce/output"
rel <- "https://github.com/yuchengwang-stat/cf-TREBLE/releases/download/v1.0.0"
sig_file  <- file.path(dir, "cfTREBLE_signature_630k_subjspec.txt.gz")
beta_file <- file.path(dir, "demo_data.beta")
options(timeout = max(3600, getOption("timeout")))
for (f in c(sig_file, beta_file))
  if (!file.exists(f)) download.file(file.path(rel, basename(f)), f, mode = "wb")

samples <- readRDS(file.path(dir, "sample_list.RDS"))
cols <- split(seq_len(sum(lengths(samples))), rep(seq_along(samples), lengths(samples)))
names(cols) <- names(samples)

signature <- as.data.frame(fread(sig_file, header = TRUE, sep = "\t"))
n <- file.size(beta_file) / 2
beta <- matrix(readBin(beta_file, "integer", 2 * n, size = 1, signed = FALSE), n, 2, byrow = TRUE)
beta <- beta[signature$index, ]

mat <- as.matrix(signature[, 4 + seq_len(sum(lengths(samples)))])
signature_list <- lapply(cols, function(idx) mat[, idx, drop = FALSE])
result <- cf_BLEND(beta[, 1], beta[, 2], signature_list)

dir.create(out, showWarnings = FALSE)
write.csv(result$cellular_frac, file.path(out, "demo_result.csv"))
written  <- read.csv(file.path(out, "demo_result.csv"), row.names = 1)
expected <- read.csv("reproduce/demo_result.csv", row.names = 1)
cat(sprintf("max absolute difference from reproduce/demo_result.csv: %.3g\n",
            max(abs(written$x - expected$x))))
print(round(result$cellular_frac[result$cellular_frac > 1e-6], 4))
