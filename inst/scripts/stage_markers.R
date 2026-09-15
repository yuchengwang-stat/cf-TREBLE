#!/usr/bin/env Rscript
## Stage 4 -- marker CpG selection for one chunk.  Array-friendly.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     chunk = Sys.getenv("SLURM_ARRAY_TASK_ID", "1"),
                     max_rows = "0", overwrite = FALSE))
cfg   <- load_config(a$config)
chunk <- as.integer(a$chunk)
tree  <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
fit   <- readRDS(file.path(run_dir(cfg, "prior", FALSE), "tree_pi.rds"))
lik   <- run_dir(cfg, "likelihood", FALSE)
dir   <- run_dir(cfg, "markers")
out   <- file.path(dir, sprintf("markers_c%03d.rds", chunk))
if (file.exists(out) && !isTRUE(a$overwrite)) { cat("already done: ", out, "\n"); quit(save = "no") }

off <- chunk_offsets(lik, as.integer(unlist(cfg$data$chunks)),
                     as.integer(cfg$data$chunk_stride %||% 900000L))
d   <- load_chunk(lik, chunk, off[[as.character(chunk)]],
                  fields = c("ll", "ll.truncate", "mu", "sigma"),
                  max_rows = as.integer(a$max_rows))
log_msg("chunk ", chunk, ": ", length(d$cpg_index), " CpGs")

## mean read depth per CpG, from the same rows of the count matrix
N <- readRDS(sprintf(cfg$data$n_chunks, chunk))[d$local_rows, , drop = FALSE]
Nmean <- rowMeans(N); rm(N); gc()

ok <- stats::complete.cases(d$fields$ll) & stats::complete.cases(d$fields$`ll.truncate`)
log_msg("usable CpGs: ", sum(ok), " (dropped ", sum(!ok), ")")
sig <- select_markers(tree, cfg$markers,
                      d$fields$ll[ok, , drop = FALSE], d$fields$`ll.truncate`[ok, , drop = FALSE],
                      d$fields$mu[ok, , drop = FALSE],  d$fields$sigma[ok, , drop = FALSE],
                      fit$pi.j, fit$pi.t, Nmean[ok], d$cpg_index[ok])
saveRDS(list(chunk = chunk, markers = sig, n_scanned = sum(ok)), out)
cat(sprintf("chunk %d: %d markers -> %s\n", chunk, if (is.null(sig)) 0L else nrow(sig), out))
