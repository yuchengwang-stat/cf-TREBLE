## Loaded by every CLI script in inst/scripts.
##
## Prefers the installed package; falls back to sourcing R/ directly so the
## scripts also run from a git checkout without installing anything.
if (requireNamespace("cfTREBLE", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cfTREBLE))
  CFT_HOME <- Sys.getenv("CFT_HOME", unset = system.file(package = "cfTREBLE"))
} else {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f[1]), "..", "..")) else getwd()
  CFT_HOME <- Sys.getenv("CFT_HOME", unset = root)
  ## Exactly the packages DESCRIPTION declares.  The original scripts also
  ## attached Matrix, expm, matrixcalc, numDeriv and fastGHQuad; none of them is
  ## called here, and requiring them would make this path fail for anyone who
  ## installed only what the README asks for.
  suppressPackageStartupMessages({
    library(data.table); library(parallel); library(pbv); library(pracma)
    library(mixsqp); library(RhpcBLASctl); library(yaml); library(jsonlite)
  })
  for (f in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)) source(f)
  message("cfTREBLE: running from the source tree at ", root)
}
## SLURM does not export OMP_NUM_THREADS.  Without it OpenMP and the BLAS see
## every core on the node while the cgroup allows only the few requested.
Sys.setenv(OMP_NUM_THREADS = "1"); RhpcBLASctl::blas_set_num_threads(1)
