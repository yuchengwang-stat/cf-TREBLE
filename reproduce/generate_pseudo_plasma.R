#!/usr/bin/env Rscript
## Generate the fragment files for one set of simulated cfDNA mixtures: a random
## set of ground-truth proportions, 20 replicates, as in the paper's simulations.
##
##   Rscript generate_pseudo_plasma.R <test_samples.rds> <pat_dir> <out_dir> <param> <nfrag> [seed]
##
## test_samples.rds  list of length 48, one entry per cell type in tree order,
##                   holding the name of the held-out test sample or NULL
## pat_dir           directory holding <sample>.hg38.pat.gz for those samples
## out_dir           where fragments_<param>_<replicate>_<j>.txt are written
## param             label for this set of proportions
## nfrag             number of fragments per mixture; sets the read depth
## seed              optional random seed
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: generate_pseudo_plasma.R <test_samples.rds> <pat_dir> <out_dir> <param> <nfrag> [seed]")
tree    <- readRDS(args[1])
pat_dir <- args[2]
out_dir <- args[3]
param   <- as.numeric(args[4])
nfrag   <- as.numeric(args[5])
if (length(args) >= 6) set.seed(as.integer(args[6]))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
idx_set <- which(!unlist(lapply(tree, is.null)))

prop <- rep(0, 48)
prop[27] <- 0.08
prop[28] <- 0.065
prop[29] <- 0.07
prop[30] <- 0.15
prop[31] <- 0.13
prop[32] <- 0.16
randomsample <- sample(idx_set[-c(22:27)], 8)
prop[randomsample[1]] <- 0.06
prop[randomsample[2]] <- 0.11
prop[randomsample[3]] <- 0.025
prop[randomsample[4]] <- 0.055
prop[randomsample[5]] <- 0.12
prop[randomsample[6]] <- 0.02
prop[randomsample[7]] <- 0.0055
prop[randomsample[8]] <- 0.0037
prop <- prop / sum(prop)
prop <- as.data.frame(t(prop))
colnames(prop) <- names(tree)
saveRDS(prop, file.path(out_dir, paste0("true_prop_", param, ".RDS")))
numbers <- floor(nfrag * prop)

j <- 1
for (i in 1:48) {
  if (is.null(tree[[i]]) || prop[i] == 0) next
  pat_path <- file.path(pat_dir, sprintf("%s.hg38.pat.gz", tree[[i]][1]))
  data <- fread(cmd = sprintf("gzip -dc %s", shQuote(pat_path)),
                sep = "\t", header = FALSE, showProgress = FALSE)
  data[, V4 := as.integer(V4)]
  data <- data[rep(seq_len(.N), data$V4)]
  data[, V4 := 1L][]
  for (k in 1:20) {
    idx <- sample(nrow(data), as.numeric(numbers[i]))
    tmp <- data[idx]
    tmp[, 4] <- 1
    fwrite(tmp, file = file.path(out_dir, sprintf("fragments_%s_%d_%d.txt", param, k, j)),
           sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
  j <- j + 1
}
