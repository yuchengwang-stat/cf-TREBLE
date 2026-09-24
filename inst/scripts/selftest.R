#!/usr/bin/env Rscript
## Self-test: runs tree -> likelihood -> prior -> markers -> signature -> deconvolve
## end to end on a small synthetic reference.  No cluster, no real data, no
## network.  Its job is to prove the wiring, not the science.
##
##   Rscript bin/selftest.R [--cpgs 120] [--tree-lists DIR] [--keep]
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(cpgs = "120", tree_lists = "", keep = FALSE, seed = "1"))
set.seed(as.integer(a$seed))
nCpG <- as.integer(a$cpgs)

tmp <- file.path(tempdir(), paste0("cft_selftest_", Sys.getpid()))
dir.create(file.path(tmp, "data"), recursive = TRUE, showWarnings = FALSE)
if (!isTRUE(a$keep)) on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
cat("workspace: ", tmp, "\n")

## ---- 1. a small tree -----------------------------------------------------------
if (nzchar(a$tree_lists)) {
  ch <- file.path(a$tree_lists, "binary.tree2.0.RDS")
  pa <- file.path(a$tree_lists, "ancestor.tree2.0.RDS")
  ct <- file.path(a$tree_lists, "all_celltype_sample_2.0.RDS")
  sl <- readRDS(ct)
} else {
  ## 8 leaves, balanced; 2 samples each
  n_leaf <- 8L
  hc <- hclust(dist(matrix(rnorm(n_leaf * 20), n_leaf)), "ward.D2")
  ctn <- paste0("CT", seq_len(n_leaf))
  hc$labels <- ctn
  tr0 <- hclust_to_tree(hc, ctn)
  sl <- setNames(lapply(ctn, function(x) paste0(x, "_", 1:2)), ctn)
  ch <- file.path(tmp, "children.rds"); saveRDS(tr0$children, ch)
  pa <- file.path(tmp, "parent.rds");   saveRDS(tr0$parent, pa)
  ct <- file.path(tmp, "celltypes.rds"); saveRDS(sl, ct)
}
saveRDS(sl, file.path(tmp, "sample_list.rds"))
groups <- sample_groups(sl); nS <- length(groups)

## ---- 2. synthetic counts with real structure -----------------------------------
## Half the CpGs are noise; the rest are markers for one randomly chosen cell type.
p_bg <- runif(nCpG, 0.2, 0.8)
P <- matrix(p_bg, nCpG, nS)
marker_of <- rep(0L, nCpG)
for (g in seq(1, nCpG, by = 2)) {
  k <- sample(seq_along(sl), 1); marker_of[g] <- k
  P[g, ] <- 0.85; P[g, groups == k] <- 0.05
}
N <- matrix(rpois(nCpG * nS, 40) + 5L, nCpG, nS)
M <- matrix(rbinom(nCpG * nS, N, P), nCpG, nS)
saveRDS(M, file.path(tmp, "data", "M_1.rds")); saveRDS(N, file.path(tmp, "data", "N_1.rds"))

## ---- 3. priors ------------------------------------------------------------------
saveRDS(list(mu = c(-1.4, 1.55, 0.66), sigma = c(0.24, 0.29, 0.65),
             lambda = c(0.06, 0.55, 0.39)), file.path(tmp, "mu_prior.rds"))
saveRDS(list(grid = c(0, .1, .15, .2, .25, .35, .5, .6, 1),
             pi = c(.052, .094, .254, .076, .251, .150, .054, .027, .041)),
        file.path(tmp, "sigma_prior.rds"))

## ---- 4. config ------------------------------------------------------------------
cfgfile <- file.path(tmp, "selftest.yaml")
writeLines(c(
  paste0("base: ", tmp), "work: work", "run_tag: selftest",
  "tree:", "  from_lists:",
  paste0("    children: ", ch), paste0("    parent: ", pa), paste0("    celltypes: ", ct),
  "priors:", paste0("  mu: ", file.path(tmp, "mu_prior.rds")),
  paste0("  sigma: ", file.path(tmp, "sigma_prior.rds")),
  "  sigma_truncate: auto", "  gh_points: 5",
  "data:", paste0("  sample_list: ", file.path(tmp, "sample_list.rds")),
  paste0("  m_chunks: ", file.path(tmp, "data", "M_%d.rds")),
  paste0("  n_chunks: ", file.path(tmp, "data", "N_%d.rds")),
  "  chunks: [1]",
  "likelihood:", paste0("  batch_size: ", max(20L, nCpG %/% 2L)), "  block_size: 10",
  paste0("  ncores: ", max(1L, min(4L, parallel::detectCores() - 1L))),
  "prior_em:", "  tol: 1.0e-5", "  max_iter: 60",
  "markers:", "  lowest_p: 0.90", "  min_read_depth: 5", "  min_rel_distance: 1.64",
  "  min_abs_distance: 0.30", "  min_abs_distance_T0: 0.30",
  "  min_IQR_target: 0.15", "  min_IQR_offtarget: 0.12", "  alpha_cap: 150000",
  "  slack: 5", "  strict_celltypes: []", "  strict_rel_distance: 1.96",
  "  subject_specific: false",
  "deconvolve:", paste0("  beta_n_cpg: ", nCpG), "  drop_sample_cols: []",
  "slurm:", "  account: none", "  partition: RM-shared", "  ntasks: 1", "  time: \"01:00:00\""
), cfgfile)

## ---- 5. run the stages ----------------------------------------------------------
R <- function(stage, ..., config = cfgfile) {
  cat("\n---- ", stage, " ----\n")
  cmd <- c(cft_script(paste0("stage_", stage, ".R")),
           "--config", config, ...)
  st <- system2("Rscript", cmd)
  if (st != 0) stop("stage ", stage, " exited ", st, call. = FALSE)
}
R("tree")
R("likelihood", "--chunk", "1", "--overwrite")
R("prior", "--chunks", "1", "--max-rows", as.character(nCpG))
R("markers", "--chunk", "1", "--overwrite")
R("signature")

## ---- 6. a synthetic plasma sample, then deconvolve ------------------------------
cfg <- load_config(cfgfile)
S   <- readRDS(file.path(run_dir(cfg, "signature", FALSE), "signature.rds"))
truth <- c(0.5, 0.2, rep(0.3 / (length(sl) - 2), length(sl) - 2))
names(truth) <- names(sl)
p_mix <- as.vector(P %*% (truth[groups] / sum(truth[groups])))
depth <- rep(60L, nCpG)
meth  <- rbinom(nCpG, depth, p_mix)
bf <- file.path(tmp, "plasma.beta")
writeBin(as.integer(t(cbind(pmin(meth, 255L), pmin(depth, 255L)))), bf, size = 1L)
R("deconvolve", "--beta", bf, "--out", file.path(tmp, "frac.csv"))

## ---- 7. assertions ---------------------------------------------------------------
cat("\n======== checks ========\n")
tr  <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
fit <- readRDS(file.path(run_dir(cfg, "prior", FALSE), "tree_pi.rds"))
est <- read.csv(file.path(tmp, "frac.csv"), row.names = 1)
chk <- function(name, ok, detail = "") cat(sprintf("  %-44s %s %s\n", name,
        if (isTRUE(ok)) "PASS" else "FAIL", detail))
pass <- TRUE; A <- function(n, ok, d = "") { chk(n, ok, d); pass <<- pass && isTRUE(ok) }

A("tree validates", tryCatch({ validate_tree(tr); TRUE }, error = function(e) FALSE))
A("n_col == n_node + #complements", tr$n_col == tr$n_node + length(tr$comp_of))
A("layer partitions internal nodes",
  setequal(unlist(tr$layer), which(lengths(tr$children) == 2L)))
A("EM produced finite pi", all(is.finite(fit$pi.j)) && all(is.finite(fit$pi.t)))
A("pi.t sums to 1", abs(sum(fit$pi.t) - 1) < 1e-6, sprintf("(%.6f)", sum(fit$pi.t)))
A("signature non-empty", nrow(S$table) > 0, sprintf("(%d CpGs)", nrow(S$table)))
A("signature beta in [0,1]", all(S$beta >= 0 & S$beta <= 1, na.rm = TRUE))
A("signature has one column per cell type", ncol(S$beta) == length(sl))
A("fractions sum to 1", abs(sum(est[1, ]) - 1) < 1e-6, sprintf("(%.6f)", sum(est[1, ])))
A("fractions non-negative", all(est[1, ] >= -1e-12))
top <- names(sl)[which.max(as.numeric(est[1, ]))]
A("dominant cell type recovered", identical(top, names(which.max(truth))),
  sprintf("(got %s, truth %s)", top, names(which.max(truth))))

cat("\n", if (pass) "SELFTEST PASSED" else "SELFTEST FAILED", "\n")
if (isTRUE(a$keep)) cat("kept workspace: ", tmp, "\n")
quit(status = if (pass) 0L else 1L)
