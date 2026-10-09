if (requireNamespace("cfTREBLE", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cfTREBLE))
  CFT_HOME <- Sys.getenv("CFT_HOME", unset = system.file(package = "cfTREBLE"))
} else {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f[1]), "..", "..")) else getwd()
  CFT_HOME <- Sys.getenv("CFT_HOME", unset = root)
  suppressPackageStartupMessages({
    library(data.table); library(parallel); library(pbv); library(pracma)
    library(mixsqp); library(RhpcBLASctl); library(yaml); library(jsonlite)
  })
  for (f in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)) source(f)
  message("cfTREBLE: running from the source tree at ", root)
}
Sys.setenv(OMP_NUM_THREADS = "1"); RhpcBLASctl::blas_set_num_threads(1)
