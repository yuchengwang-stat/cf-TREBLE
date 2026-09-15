#!/usr/bin/env Rscript
## Stage 6 -- deconvolve one or more .beta samples against the signature.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     beta = "", signature = "", out = "", celltypes = "",
                     subject = "", na = "drop"))
cfg <- load_config(a$config)
if (!nzchar(a$beta)) die("--beta <file[,file...] | glob> is required")
sigfile <- if (nzchar(a$signature)) a$signature else file.path(run_dir(cfg, "signature", FALSE), "signature.rds")
S <- read_signature(sigfile, if (nzchar(a$celltypes)) strsplit(a$celltypes, ",")[[1]] else NULL)

## `drop_sample_cols` is expressed against the reference panel, so it applies
## only to a signature this pipeline produced, and only removes a cell type whose
## samples are all dropped -- the beta values themselves are already per cell type.
drop <- as.integer(unlist(cfg$deconvolve$drop_sample_cols))
if (length(drop) && grepl("[.]rds$", sigfile, ignore.case = TRUE) &&
    !is.null(cfg$data$sample_list) && file.exists(cfg$data$sample_list)) {
  sl <- readRDS(cfg$data$sample_list)
  gone <- names(sl)[!seq_along(sl) %in% unique(sample_groups(sl)[-drop])]
  gone <- intersect(gone, colnames(S$beta))
  if (length(gone)) {
    log_msg("dropping cell types with no reference sample left: ",
            paste(gone, collapse = ", "))
    S$beta <- S$beta[, setdiff(colnames(S$beta), gone), drop = FALSE]
  }
}
## Without --subject each cell type is one column, its own signature value.  With
## it, each cell type is its reference panel: one column per subject, which is the
## mixture BLEND actually models.
idx <- S$index
## --subject on its own uses this run's own stage-`subject` output; give it a
## directory to point somewhere else.
use_subject <- isTRUE(a$subject) || (is.character(a$subject) && nzchar(a$subject))
if (use_subject) {
  subdir <- if (is.character(a$subject) && dir.exists(a$subject)) a$subject
            else run_dir(cfg, "subject", FALSE)
  sl <- readRDS(cfg$data$sample_list)
  sub <- read_subject(subdir, idx, na = a$na)
  P <- sub$beta
  nas <- !stats::complete.cases(P)
  if (any(nas)) {
    if (identical(a$na, "drop")) {
      log_msg("dropping ", sum(nas), " of ", nrow(P),
              " CpGs where some reference subject has no reads")
      P <- P[!nas, , drop = FALSE]; idx <- idx[!nas]
      S$beta <- S$beta[!nas, , drop = FALSE]
    } else {
      ct <- S$beta[, sample_groups(sl), drop = FALSE]
      hole <- is.na(P); P[hole] <- ct[hole]
      log_msg("filled ", sum(hole), " subject/CpG holes with the cell type's own value")
    }
  }
  phi <- subject_phi(P, sl, as.integer(unlist(cfg$deconvolve$drop_sample_cols)))
  log_msg("subject-specific signature: ", nrow(P), " CpGs, ", length(phi),
          " cell types, ", sum(vapply(phi, ncol, 0L)), " reference samples")
} else {
  phi <- signature_phi(S$beta)
  log_msg("signature: ", nrow(S$beta), " CpGs, ", length(phi), " cell types (",
          basename(sigfile), ")")
}

files <- unlist(lapply(strsplit(a$beta, ",")[[1]], function(p)
  if (grepl("[*?]", p)) Sys.glob(p) else p))
if (!length(files)) die("no .beta files matched ", a$beta)

res <- t(vapply(files, function(f) {
  log_msg("deconvolving ", basename(f))
  deconvolve_beta(f, if (is.null(cfg$deconvolve$beta_n_cpg)) NULL else
                     as.integer(cfg$deconvolve$beta_n_cpg), idx, phi)
}, numeric(length(phi))))
rownames(res) <- basename(files)

out <- if (nzchar(a$out)) a$out else file.path(run_dir(cfg, "deconvolve"), "cell_fractions.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out)
write_manifest(run_dir(cfg, "deconvolve"), "deconvolve",
               entries = list(list(path = out)),
               extra = list(samples = rownames(res), n_cpg = length(idx),
                            celltypes = colnames(res), signature = sigfile,
                            subject = use_subject, n_ref = sum(vapply(phi, ncol, 0L))))
cat("wrote ", out, "\n"); print(round(res, 4))
