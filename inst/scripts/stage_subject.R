#!/usr/bin/env Rscript
## Stage 7 -- subject-specific signature for one chunk's marker CpGs.
##
## Runs after `signature`.  Array-friendly over chunks.  The 19,458-wide
## per-CpG intermediate is consumed in the same pass and never written: the
## output is markers x samples.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     chunk = Sys.getenv("SLURM_ARRAY_TASK_ID", "1"),
                     coverage = "", exact_sample = "0", max_cpgs = "0", overwrite = FALSE))
cfg   <- load_config(a$config)
chunk <- as.integer(a$chunk)
cov   <- if (nzchar(a$coverage)) as.numeric(a$coverage) else
         as.numeric(cfg$subject$weight_coverage %||% 0.9999)
tree  <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
fit   <- readRDS(file.path(run_dir(cfg, "prior", FALSE), "tree_pi.rds"))
lik   <- run_dir(cfg, "likelihood", FALSE)
dir   <- run_dir(cfg, "subject")
out   <- file.path(dir, sprintf("subject_c%03d.rds", chunk))
if (file.exists(out) && !isTRUE(a$overwrite)) { cat("already done: ", out, "\n"); quit(save = "no") }

mfile <- file.path(run_dir(cfg, "markers", FALSE), sprintf("markers_c%03d.rds", chunk))
mk <- readRDS(mfile)$markers
if (is.null(mk) || !nrow(mk)) { cat("no markers in chunk ", chunk, "\n"); quit(save = "no") }

off <- chunk_offsets(lik, as.integer(unlist(cfg$data$chunks)),
                     as.integer(cfg$data$chunk_stride %||% 900000L))
d <- load_chunk(lik, chunk, off[[as.character(chunk)]],
                fields = c("ll", "ll.truncate"))
rows <- match(mk$index, d$cpg_index)
keep <- !is.na(rows); rows <- rows[keep]; mk <- mk[keep, , drop = FALSE]
if (as.integer(a$max_cpgs) > 0L && as.integer(a$max_cpgs) < length(rows)) {
  sel <- seq_len(as.integer(a$max_cpgs)); rows <- rows[sel]; mk <- mk[sel, , drop = FALSE]
}
log_msg("chunk ", chunk, ": ", length(rows), " marker CpGs, weight coverage ", cov)

## p2 / p.t at the marker rows only
pa  <- calculate.p.gj(tree$children, fit$pi.j, d$fields$ll[rows, , drop = FALSE], tree$layer)
p.t <- calculate.p.t(fit$pi.t, fit$pi.j, tree, d$fields$ll[rows, , drop = FALSE],
                     d$fields$`ll.truncate`[rows, , drop = FALSE], pa[[2]])
p2  <- node_joint_p(tree, pa[[1]]) * p.t[, 1]
rm(pa, d); gc()

sl <- readRDS(cfg$data$sample_list); groups <- sample_groups(sl)
pr <- load_priors(cfg)
plan <- subject_weight_plan(tree)
loc <- mk$local_idx
M <- readRDS(sprintf(cfg$data$m_chunks, chunk))[loc, , drop = FALSE]
N <- readRDS(sprintf(cfg$data$n_chunks, chunk))[loc, , drop = FALSE]
U <- N - M; rm(N); gc()

ncores <- as.integer(cfg$likelihood$ncores %||% 1L)
run <- function(idx, coverage) {
  blocks <- split(idx, ceiling(seq_along(idx) / max(1L, ceiling(length(idx) / ncores))))
  res <- parallel::mclapply(blocks, function(b) lapply(b, function(g)
      subject_one_cpg(tree, plan, groups, M[g, ], U[g, ], p2[g, ], p.t[g, ], pr, coverage)),
    mc.cores = ncores, mc.preschedule = FALSE)
  unlist(res, recursive = FALSE)
}

t0 <- Sys.time()
r <- run(seq_along(rows), cov)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
P <- do.call(rbind, lapply(r, `[[`, "p"))
colnames(P) <- unlist(sl, use.names = FALSE)
nodes_used <- vapply(r, `[[`, 0L, "n_nodes")
cover_real <- do.call(rbind, lapply(r, `[[`, "covered"))
## A zero here means the sample had no reads at that CpG, not a pruning failure,
## so report the minimum over the sample/CpG pairs that carry any weight at all.
nz <- cover_real[cover_real > 1e-12]
log_msg(sprintf("pruned pass: %.1f min, nodes/CpG median %d of %d, realised coverage min %.6f over %d of %d sample-CpG pairs with data",
                el, median(nodes_used), tree$n_col,
                if (length(nz)) min(nz) else NA_real_, length(nz), length(cover_real)))

## Optional: redo a sample of CpGs with every node, and report the real error.
acc <- NULL
ns <- as.integer(a$exact_sample)
if (ns > 0L && cov < 1) {
  ns <- min(ns, length(rows))
  set.seed(1); idx <- sort(sample(seq_along(rows), ns))
  t1 <- Sys.time(); re <- run(idx, 1.0)
  el2 <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
  E <- do.call(rbind, lapply(re, `[[`, "p"))
  dif <- abs(P[idx, , drop = FALSE] - E)
  dif[is.na(dif)] <- 0   # no-data samples are NA in both and carry no error
  acc <- list(n = ns, max = max(dif, na.rm = TRUE), mean = mean(dif, na.rm = TRUE),
              q999 = quantile(dif, .999, na.rm = TRUE), bound = 1 - cov,
              min_per_cpg = min(cover_real[idx, ], na.rm = TRUE),
              exact_min = el2, pruned_min = el * ns / length(rows))
  log_msg(sprintf("exact check on %d CpGs: max |diff| %.3e (bound %.1e), mean %.3e; exact took %.1f min vs %.2f min pruned (%.1fx)",
                  ns, acc$max, acc$bound, acc$mean, el2, acc$pruned_min, el2 / max(acc$pruned_min, 1e-9)))
}

saveRDS(list(chunk = chunk, index = mk$index, celltype = mk$celltype,
             p_subj = P, samples = colnames(P), coverage = cov,
             nodes_used = nodes_used, accuracy = acc, minutes = el), out)
write_manifest(dir, "subject", entries = list(list(path = out)),
               extra = list(chunk = chunk, n_cpg = nrow(P), coverage = cov,
                            median_nodes = median(nodes_used), minutes = el,
                            accuracy = acc))
cat(sprintf("chunk %d: %d CpGs x %d subjects -> %s\n", chunk, nrow(P), ncol(P), out))
