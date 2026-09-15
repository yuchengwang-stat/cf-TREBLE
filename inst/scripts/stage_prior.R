#!/usr/bin/env Rscript
## Stage 3 -- EM for the Bernoulli / categorical priors on the tree.
##
## Runs over EVERY CpG, not a sample.  pi.t has one component per partition and
## the rare ones are what marker selection keys on, so their proportions have to
## be estimated from the whole genome; a subsample would be noisy exactly where
## it matters.  The EM update is a sum over CpGs, so chunks are accumulated one
## at a time and the 83 GB of likelihood never has to be resident at once --
## bin/test_em_streaming.R asserts this is identical to the batch computation,
## not an approximation.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     chunks = "", cache = "", max_rows = "0"))
cfg  <- load_config(a$config)
tree <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
lik  <- run_dir(cfg, "likelihood", FALSE)
dir  <- run_dir(cfg, "prior")

chunks <- if (nzchar(a$chunks)) as.integer(strsplit(a$chunks, ",")[[1]]) else
          as.integer(unlist(cfg$data$chunks))
have <- chunks[vapply(chunks, function(k)
  length(list.files(lik, pattern = sprintf("^lik_c%03d_b", k))) > 0, TRUE)]
if (!length(have)) die("no likelihood output found in ", lik)
if (length(have) < length(chunks))
  log_msg("WARNING: only ", length(have), " of ", length(chunks),
          " chunks are on disk -- the prior will describe those CpGs only")
log_msg("EM over chunks: ", paste(have, collapse = ","))

off <- chunk_offsets(lik, as.integer(unlist(cfg$data$chunks)),
                     as.integer(cfg$data$chunk_stride %||% 900000L))
max_rows <- as.integer(a$max_rows)
loader <- function(k) {
  cd <- load_chunk(lik, k, off[[as.character(k)]], fields = c("ll", "ll.truncate"),
                   max_rows = max_rows)
  ok <- stats::complete.cases(cd$fields$ll) & stats::complete.cases(cd$fields$`ll.truncate`)
  if (!any(ok)) return(NULL)
  list(ll = cd$fields$ll[ok, , drop = FALSE], llt = cd$fields$`ll.truncate`[ok, , drop = FALSE])
}
cache <- if (nzchar(a$cache)) isTRUE(as.logical(a$cache)) else isTRUE(cfg$prior_em$cache)
if (cache) log_msg("caching every chunk in memory (~2.7 GB per 900k-row chunk)")

t0 <- Sys.time()
fit <- estimate_tree_prior_streaming(tree, have, loader,
                                     tol = as.numeric(cfg$prior_em$tol),
                                     max_iter = as.integer(cfg$prior_em$max_iter),
                                     cache = cache)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
save_rds(fit, file.path(dir, "tree_pi.rds"))
write_manifest(dir, "prior", entries = list(list(path = file.path(dir, "tree_pi.rds"))),
               extra = list(chunks = have, n_cpg = fit$n_cpg, iters = length(fit$trace),
                            converged = fit$converged, final_err = tail(fit$trace, 1),
                            minutes = el, streaming = TRUE))
cat(sprintf("tree prior: %d CpGs, %d iters, final err %.3e, converged=%s, %.1f min\n",
            fit$n_cpg, length(fit$trace), tail(fit$trace, 1), fit$converged, el))
