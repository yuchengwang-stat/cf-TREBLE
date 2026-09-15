#!/usr/bin/env Rscript
## `cftreble submit` builds the SLURM chain correctly, checked without a cluster
## by putting a stub `sbatch` on PATH and reading back what it was called with.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "..", "scripts", "boot.R"))
tmp <- file.path(tempdir(), paste0("cft_submit_", Sys.getpid()))
dir.create(file.path(tmp, "stub"), recursive = TRUE, showWarnings = FALSE)

log <- file.path(tmp, "calls.log")
writeLines(c("#!/bin/bash", sprintf('echo "$*" >> %s', shQuote(log)), "echo 111222"),
           sb <- file.path(tmp, "stub", "sbatch"))
Sys.chmod(sb, "0755")

## a miniature dataset to submit against
ex <- file.path(tmp, "example")
system2("Rscript", c(cft_script("stage_example.R"), "--out", ex, "--cpgs", "200", "--chunks", "3"),
        stdout = FALSE, stderr = FALSE)
cfgfile <- file.path(ex, "example.yaml")
stopifnot(file.exists(cfgfile))

old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(file.path(tmp, "stub"), old_path, sep = .Platform$path.sep))
out <- system2(cft_script("cftreble"), c("--config", cfgfile, "submit"),
               stdout = TRUE, stderr = TRUE)
Sys.setenv(PATH = old_path)

calls <- if (file.exists(log)) readLines(log) else character(0)
stage_of <- function(x) sub(".*run_stage[.]sbatch ([a-z]+) .*", "\\1", x)
got <- vapply(calls, stage_of, "")
cat("  stages submitted:", paste(got, collapse = " -> "), "\n")
stopifnot(identical(unname(got), c("likelihood", "prior", "markers", "signature", "subject")))

## the per-chunk stages are arrays over exactly the configured chunks
arr <- grepl("--array=1-3", calls)
cat(sprintf("  array over chunks 1-3 on: %s\n", paste(got[arr], collapse = ", ")))
stopifnot(identical(unname(got[arr]), c("likelihood", "markers", "subject")))
stopifnot(!any(arr[got %in% c("prior", "signature")]))

## every stage waits for the one before it
dep <- sub(".*--dependency=afterok:([0-9]+).*", "\\1", calls)
stopifnot(dep[1] == calls[1],                       # likelihood has no dependency
          all(dep[-1] == "111222"))                 # the rest chain off the stub's id
cat("  dependency chain present on prior, markers, signature, subject   PASS\n")

## resources come from the config, not from the script
stopifnot(all(grepl("-p RM-shared", calls)), all(grepl("--ntasks-per-node=", calls)))
cat("  partition and ntasks taken from the config                       PASS\n")

## and stage 6 is skipped when the config says so
cfg2 <- file.path(tmp, "nosubject.yaml")
writeLines(sub("^  subject_specific:.*$", "  subject_specific:      false",
               readLines(cfgfile)), cfg2)
unlink(log)
Sys.setenv(PATH = paste(file.path(tmp, "stub"), old_path, sep = .Platform$path.sep))
o2 <- system2(cft_script("cftreble"), c("--config", cfg2, "submit"), stdout = TRUE, stderr = TRUE)
Sys.setenv(PATH = old_path)
got2 <- vapply(readLines(log), stage_of, "")
stopifnot(!"subject" %in% got2, any(grepl("subject    -> skipped", o2)))
cat("  subject_specific: false omits stage 6                            PASS\n")

unlink(tmp, recursive = TRUE)
cat("\nSUBMIT TESTS PASSED\n")
