## cf-TREBLE pipeline -- shared utilities -------------------------------------
suppressPackageStartupMessages({ library(yaml); library(jsonlite) })

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     paste0(...)), file = stderr())
die <- function(...) stop(paste0(...), call. = FALSE)

## ---- config -----------------------------------------------------------------
## Resolves ${key} references against already-resolved top-level scalars.
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  scal <- cfg[vapply(cfg, function(x) is.character(x) && length(x) == 1, TRUE)]
  subst <- function(x) {
    if (is.character(x)) {
      for (k in names(scal)) x <- gsub(paste0("${", k, "}"), scal[[k]], x, fixed = TRUE)
      x
    } else if (is.list(x)) lapply(x, subst) else x
  }
  cfg <- subst(cfg)
  cfg$base <- normalizePath(cfg$base, mustWork = FALSE)
  cfg$work <- if (startsWith(cfg$work, "/")) cfg$work else file.path(cfg$base, cfg$work)
  cfg$.file <- normalizePath(path, mustWork = FALSE)
  cfg
}

## run directory: work/<run_tag>/<stage>
run_dir <- function(cfg, stage, create = TRUE) {
  d <- file.path(cfg$work, cfg$run_tag, stage)
  if (create) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## ---- manifests ---------------------------------------------------------------
## A stage writes one JSON manifest naming every artefact it produced, so the
## next stage never has to guess a filename or a row order.
write_manifest <- function(dir, stage, entries, extra = list()) {
  man <- c(list(stage = stage,
                created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                host = Sys.info()[["nodename"]],
                entries = entries), extra)
  p <- file.path(dir, "manifest.json")
  write_json(man, p, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  log_msg("manifest -> ", p)
  invisible(p)
}
read_manifest <- function(dir) {
  p <- file.path(dir, "manifest.json")
  if (!file.exists(p)) die("no manifest in ", dir, " -- run the previous stage first")
  fromJSON(p, simplifyVector = FALSE)
}

## ---- misc --------------------------------------------------------------------
## sample_list -> integer vector mapping each column of M/N to its cell type
sample_groups <- function(sample_list) {
  rep(seq_along(sample_list), times = lengths(sample_list))
}

## Never silently append to an existing artefact (the old pipeline's worst trap).
fresh_file <- function(path, overwrite) {
  if (file.exists(path) && !overwrite)
    die(path, " already exists; pass --overwrite to replace it")
  if (file.exists(path)) unlink(path)
  path
}

save_rds <- function(obj, path) { saveRDS(obj, path); log_msg("wrote ", path, " (",
    format(structure(file.size(path), class = "object_size"), units = "auto"), ")"); path }

args_parse <- function(defaults = list()) {
  a <- commandArgs(trailingOnly = TRUE)
  out <- defaults; i <- 1
  while (i <= length(a)) {
    k <- a[i]
    if (!startsWith(k, "--")) die("unexpected argument: ", k)
    k <- sub("^--", "", k)
    if (i < length(a) && !startsWith(a[i + 1], "--")) { out[[gsub("-", "_", k)]] <- a[i + 1]; i <- i + 2 }
    else { out[[gsub("-", "_", k)]] <- TRUE; i <- i + 1 }
  }
  out
}

## ---- where the packaged files live ---------------------------------------------
## Work from an installed package or from a git checkout, without the caller
## having to know which.
cft_file <- function(...) {
  p <- system.file(..., package = "cfTREBLE")
  if (nzchar(p)) return(p)
  root <- Sys.getenv("CFT_HOME", unset = getwd())
  file.path(root, "inst", ...)
}
cft_default_config <- function() cft_file("config", "default.yaml")
cft_script <- function(name) cft_file("scripts", name)
cft_slurm  <- function(name) cft_file("slurm", name)
