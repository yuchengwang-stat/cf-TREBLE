#!/usr/bin/env Rscript
## Validate everything a config points at, before spending compute on it.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(), chunks = "", deep = FALSE))
cfg <- load_config(a$config)
cat("checking ", a$config, "\n", sep = "")
ch <- if (nzchar(a$chunks)) as.integer(strsplit(a$chunks, ",")[[1]]) else NULL
res <- check_inputs(cfg, chunks = ch, deep = isTRUE(a$deep))
quit(status = if (any(res$result != "ok")) 1L else 0L)
