#!/usr/bin/env Rscript
## Stage 2 -- per-node log-likelihood for one CpG chunk.  Array-friendly:
## --chunk comes from SLURM_ARRAY_TASK_ID.  Finished batches are skipped.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     chunk = Sys.getenv("SLURM_ARRAY_TASK_ID", "1"),
                     fields = "min", max_rows = "0", sample_rows = "0",
                     sample_seed = "42", overwrite = FALSE))
cfg   <- load_config(a$config)
chunk <- as.integer(a$chunk)
dir   <- run_dir(cfg, "likelihood")
tree  <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
pr    <- load_priors(cfg)
fields <- switch(a$fields, min = LIK_FIELDS_MIN, all = LIK_FIELDS_ALL,
                 strsplit(a$fields, ",")[[1]])
sl    <- readRDS(cfg$data$sample_list)
groups <- sample_groups(sl)

Mf <- sprintf(cfg$data$m_chunks, chunk); Nf <- sprintf(cfg$data$n_chunks, chunk)
log_msg("chunk ", chunk, ": ", basename(Mf))
M <- readRDS(Mf); N <- readRDS(Nf)
if (ncol(M) != length(groups))
  die("M has ", ncol(M), " columns, sample_list implies ", length(groups))
## Record the FULL height of the chunk before any truncation: global CpG ids of
## later chunks are derived from it, and a --max-rows run must not shift them.
writeLines(as.character(nrow(M)), file.path(dir, sprintf("rows_c%03d.txt", chunk)))
## --max-rows takes the HEAD of the chunk, which is one contiguous stretch of one
## chromosome: fine for a wiring check, useless for judging a signature, because
## cell-type markers are scattered genome-wide.  --sample-rows draws at random
## instead, so a small test still sees the whole chunk.
sub_rows <- NULL
max_rows <- as.integer(a$max_rows); samp_rows <- as.integer(a$sample_rows)
if (samp_rows > 0L && samp_rows < nrow(M)) {
  set.seed(as.integer(a$sample_seed) + chunk)
  sub_rows <- sort(sample.int(nrow(M), samp_rows))
  log_msg("sampling ", samp_rows, " rows at random from ", nrow(M))
} else if (max_rows > 0L && max_rows < nrow(M)) {
  sub_rows <- seq_len(max_rows)
  log_msg("taking the first ", max_rows, " rows (contiguous -- wiring checks only)")
}
if (!is.null(sub_rows)) { M <- M[sub_rows, , drop = FALSE]; N <- N[sub_rows, , drop = FALSE] }

## ---- pre-filter -------------------------------------------------------------
## The gate is marker selection's own requirement, read from its config key, so
## the two cannot drift apart.  Rows that fail it are skipped -- except for a
## uniform random sample, which stage `prior` uses so the EM stays unbiased.
gate_rows <- seq_len(nrow(M)); em_rows <- integer(0)
if (isTRUE(cfg$likelihood$prefilter)) {
  thr <- as.numeric(cfg$markers$min_read_depth)
  Nmean <- rowMeans(N)
  gate_rows <- which(Nmean >= thr)
  ns <- min(as.integer(cfg$likelihood$em_sample_per_chunk %||% 0L), nrow(M))
  if (ns > 0L) { set.seed(1000L + chunk); em_rows <- sort(sample.int(nrow(M), ns)) }
  keep <- sort(union(gate_rows, em_rows))
  log_msg(sprintf("prefilter mean depth >= %g: %d of %d rows (%.1f%%), plus %d EM rows -> %d computed",
                  thr, length(gate_rows), nrow(M), 100 * length(gate_rows) / nrow(M),
                  length(setdiff(em_rows, gate_rows)), length(keep)))
  saveRDS(list(chunk = chunk, n_total = nrow(M), gate_rows = gate_rows,
               em_rows = em_rows, threshold = thr),
          file.path(dir, sprintf("prefilter_c%03d.rds", chunk)))
  M <- M[keep, , drop = FALSE]; N <- N[keep, , drop = FALSE]
  attr(M, "orig_rows") <- keep
}
gate_keep <- attr(M, "orig_rows")
## Map back to positions in the FULL chunk, so cpg_index stays a .beta row number
## whatever subsetting happened above.
orig_rows <- if (is.null(sub_rows)) {
  if (is.null(gate_keep)) seq_len(nrow(M)) else gate_keep
} else {
  if (is.null(gate_keep)) sub_rows else sub_rows[gate_keep]
}

bs <- as.integer(cfg$likelihood$batch_size)
nb <- ceiling(nrow(M) / bs)

for (b in seq_len(nb)) {
  out <- lik_path(dir, chunk, b)
  rng <- seq.int((b - 1L) * bs + 1L, min(b * bs, nrow(M)))
  ## Resuming must check WHAT is in the file, not just that it exists.  A batch
  ## written by an earlier --max-rows run covers fewer rows than this one wants;
  ## reusing it would silently drop CpGs from the signature.
  if (file.exists(out) && !isTRUE(a$overwrite)) {
    prev <- tryCatch(readRDS(out), error = function(e) NULL)
    if (!is.null(prev) && identical(as.integer(prev$local_rows), as.integer(orig_rows[rng])) &&
        identical(as.integer(prev$n_col), as.integer(tree$n_col)) &&
        all(fields %in% names(prev$fields))) {
      log_msg("batch ", b, "/", nb, " already covers the same rows, skip")
      next
    }
    log_msg("batch ", b, "/", nb, " on disk does not match this run (",
            if (is.null(prev)) "unreadable" else
              paste0("has ", length(prev$local_rows), " rows, want ", length(rng)),
            ") -- recomputing")
  }
  t0 <- Sys.time()
  res <- lik_batch(M[rng, , drop = FALSE], N[rng, , drop = FALSE], tree, groups, pr,
                   fields, as.integer(cfg$likelihood$ncores), as.integer(cfg$likelihood$block_size))
  saveRDS(list(chunk = chunk, batch = b, local_rows = orig_rows[rng],
               fields = res, n_col = tree$n_col), out)
  log_msg(sprintf("batch %d/%d  rows=%d  %.1f min -> %s",
                  b, nb, length(rng), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  basename(out)))
}
done <- list.files(dir, pattern = sprintf("^lik_c%03d_b.*rds$", chunk), full.names = TRUE)
cat(sprintf("chunk %d complete: %d/%d batches\n", chunk, length(done), nb))
