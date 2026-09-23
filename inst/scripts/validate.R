#!/usr/bin/env Rscript
## Accuracy check: build a reference whose truth we know, run the whole pipeline
## on it, then deconvolve plasma samples mixed at known proportions and measure
## how close the estimates are.
##
##   Rscript bin/validate.R [--cpgs 3000] [--depth 40] [--plasma-depth 30]
##                          [--tree-lists DIR] [--seed 1] [--keep]
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(cpgs = "3000", depth = "40", plasma_depth = "30",
                     tree_lists = "", seed = "1", keep = FALSE, ncores = ""))
set.seed(as.integer(a$seed))
nCpG <- as.integer(a$cpgs)
ncores <- if (nzchar(a$ncores)) as.integer(a$ncores) else
  max(1L, min(26L, parallel::detectCores() - 1L))

tmp <- file.path(Sys.getenv("CFT_SCRATCH", tempdir()),
                 paste0("cft_validate_", Sys.getpid()))
dir.create(file.path(tmp, "data"), recursive = TRUE, showWarnings = FALSE)
if (!isTRUE(a$keep)) on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
cat("workspace: ", tmp, "\n  CpGs=", nCpG, " ref depth=", a$depth,
    " plasma depth=", a$plasma_depth, " cores=", ncores, "\n", sep = "")

## ---- tree + sample layout --------------------------------------------------------
if (nzchar(a$tree_lists)) {
  ch <- file.path(a$tree_lists, "binary.tree2.0.RDS")
  pa <- file.path(a$tree_lists, "ancestor.tree2.0.RDS")
  ctf <- file.path(a$tree_lists, "all_celltype_sample_2.0.RDS")
  sl <- readRDS(ctf)
  tree <- tree_from_lists(ch, pa, ctf)
} else {
  n_leaf <- 16L
  hc <- hclust(dist(matrix(rnorm(n_leaf * 30), n_leaf)), "ward.D2")
  ctn <- sprintf("CT%02d", seq_len(n_leaf)); hc$labels <- ctn
  tree <- hclust_to_tree(hc, ctn)
  sl <- setNames(lapply(ctn, function(x) paste0(x, "_", 1:3)), ctn)
  ch <- file.path(tmp, "children.rds"); saveRDS(tree$children, ch)
  pa <- file.path(tmp, "parent.rds");   saveRDS(tree$parent, pa)
  ctf <- file.path(tmp, "celltypes.rds"); saveRDS(sl, ctf)
}
saveRDS(sl, file.path(tmp, "sample_list.rds"))
groups <- sample_groups(sl)
cat("tree: ", tree$n_leaf, " cell types, ", tree$n_node, " nodes, ",
    length(groups), " reference samples\n", sep = "")

## ---- truth, reference, plasma -----------------------------------------------------
truth <- simulate_truth(tree, nCpG)
cat("truth: ", sum(truth$kind == "onevrest"), " one-vs-rest CpGs, ",
    sum(truth$kind == "class"), " class CpGs, ",
    sum(truth$kind == "bg"), " background\n", sep = "")
ref <- simulate_reference(truth, groups, depth = as.integer(a$depth))
saveRDS(ref$M, file.path(tmp, "data", "M_1.rds"))
saveRDS(ref$N, file.path(tmp, "data", "N_1.rds"))

saveRDS(list(mu = c(-1.4216, 1.5526, 0.6610), sigma = c(0.2365, 0.2875, 0.6497),
             lambda = c(0.0580, 0.5485, 0.3935)), file.path(tmp, "mu_prior.rds"))
saveRDS(list(grid = c(0, .10, .15, .20, .25, .35, .50, .60, 1.00),
             pi = c(.0522, .0943, .2542, .0764, .2508, .1504, .0544, .0268, .0405)),
        file.path(tmp, "sigma_prior.rds"))

cfgfile <- file.path(tmp, "validate.yaml")
writeLines(c(
  paste0("base: ", tmp), "work: work", "run_tag: validate",
  "tree:", "  from_lists:", paste0("    children: ", ch),
  paste0("    parent: ", pa), paste0("    celltypes: ", ctf),
  "priors:", paste0("  mu: ", file.path(tmp, "mu_prior.rds")),
  paste0("  sigma: ", file.path(tmp, "sigma_prior.rds")),
  "  sigma_truncate: auto", "  gh_points: 5",
  "data:", paste0("  sample_list: ", file.path(tmp, "sample_list.rds")),
  paste0("  m_chunks: ", file.path(tmp, "data", "M_%d.rds")),
  paste0("  n_chunks: ", file.path(tmp, "data", "N_%d.rds")),
  "  chunks: [1]", paste0("  chunk_stride: ", nCpG),
  "likelihood:", paste0("  batch_size: ", nCpG), "  block_size: 2000",
  paste0("  ncores: ", ncores),
  "prior_em:", "  tol: 1.0e-6", "  max_iter: 400",
  "markers:", "  lowest_p: 0.90", "  min_read_depth: 15", "  min_rel_distance: 1.64",
  "  min_abs_distance: 0.30", "  min_abs_distance_T0: 0.30",
  "  min_IQR_target: 0.15", "  min_IQR_offtarget: 0.12", "  alpha_cap: 150000",
  "  slack: 5", "  strict_celltypes: []", "  strict_rel_distance: 1.96",
  "  strict_max_samples: 0", "  subject_specific: false",
  "deconvolve:", "  beta_n_cpg: null", "  drop_sample_cols: []",
  "slurm:", "  account: none", "  partition: RM-shared", "  ntasks: 1", "  time: \"01:00:00\""
), cfgfile)

R <- function(stage, ...) {
  st <- system2("Rscript", c(cft_script(paste0("stage_", stage, ".R")),
                             "--config", cfgfile, ...), stdout = NULL, stderr = NULL)
  if (st != 0) stop("stage ", stage, " exited ", st, call. = FALSE)
  cat("  [", stage, " ok]\n", sep = "")
}
cat("\nrunning pipeline ...\n")
t0 <- Sys.time()
R("tree"); R("likelihood", "--chunk", "1", "--overwrite")
R("prior", "--chunks", "1", "--max-rows", as.character(nCpG))
R("markers", "--chunk", "1", "--overwrite"); R("signature")
cat(sprintf("  pipeline wall time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cfg <- load_config(cfgfile)
S <- readRDS(file.path(run_dir(cfg, "signature", FALSE), "signature.rds"))

## ---- did marker selection find the right CpGs? ------------------------------------
cat("\n======== marker selection vs truth ========\n")
sel <- S$table$index
is_marker <- truth$kind != "bg"
tp <- sum(is_marker[sel]); fp <- sum(!is_marker[sel])
cat(sprintf("  selected            %d of %d CpGs\n", length(sel), nCpG))
cat(sprintf("  precision           %.3f   (%d true markers, %d background)\n",
            tp / max(1, length(sel)), tp, fp))
cat(sprintf("  recall              %.3f   (of %d planted markers)\n",
            tp / max(1, sum(is_marker)), sum(is_marker)))
## The signature table lists targets by label, and drops the node index; the
## per-chunk marker table keeps it, so score against that.  A CpG can be selected
## for several nodes, and it counts as correct if the planted one is among them.
mk <- readRDS(file.path(run_dir(cfg, "markers", FALSE), "markers_c001.rds"))
mk <- if (is.data.frame(mk)) mk else mk$markers
nodes_of_cpg <- split(as.integer(mk$node), mk$index)
true_sel <- sel[is_marker[sel]]
right_node <- sum(vapply(as.character(true_sel), function(i)
  !is.null(nodes_of_cpg[[i]]) && truth$node_of[as.integer(i)] %in% nodes_of_cpg[[i]], TRUE))
cat(sprintf("  correct node        %.3f   (of the %d true markers selected)\n",
            right_node / max(1, tp), tp))

## ---- how close are the estimated fractions? ---------------------------------------
cat("\n======== deconvolution vs known proportions ========\n")
phi <- signature_to_phi(S$beta[, groups, drop = FALSE], sl, integer(0))
mix <- make_mixtures(tree$celltypes)
rows <- list()
for (nm in names(mix)) {
  bf <- simulate_plasma(truth, mix[[nm]], as.integer(a$plasma_depth),
                        file.path(tmp, paste0("plasma_", gsub("[^A-Za-z0-9]", "_", nm), ".beta")))
  est <- deconvolve_beta(bf, NULL, sel, phi)
  m <- deconv_metrics(est, mix[[nm]])
  rows[[nm]] <- m
  cat(sprintf("\n  %-16s RMSE %.4f  MAE %.4f  r %.3f  max err %.3f\n",
              nm, m["RMSE"], m["MAE"], m["pearson"], m["max_abs_err"]))
  ## Only meaningful when the truth actually has a unique largest component.
  tied <- sum(mix[[nm]] == max(mix[[nm]])) > 1L
  if (tied) {
    cat(sprintf("  %-16s dominant: truth is a %d-way tie, skipped\n", "",
                sum(mix[[nm]] == max(mix[[nm]]))))
  } else {
    top_t <- tree$celltypes[which.max(mix[[nm]])]; top_e <- names(which.max(est))
    cat(sprintf("  %-16s dominant: truth %s / est %s %s\n", "",
                top_t, top_e, if (identical(top_t, top_e)) "OK" else "MISMATCH"))
  }
  ## For a spiked mixture the number that matters is the spiked type itself.
  sp <- attr(mix[[nm]], "spike")
  if (!is.null(sp)) {
    sf <- attr(mix[[nm]], "spike_frac"); se <- est[[sp]]
    cat(sprintf("  %-16s SPIKE %s: truth %.4f  est %.4f  (%s; background max %.4f)\n", "",
                sp, sf, se,
                if (se >= sf / 3) "detected" else "MISSED",
                max(est[setdiff(names(est), sp)])))
  }
  ord <- order(mix[[nm]], decreasing = TRUE)[1:min(5, length(est))]
  for (i in ord) if (mix[[nm]][i] > 0 || est[i] > 0.01)
    cat(sprintf("       %-22s truth %.4f  est %.4f\n", tree$celltypes[i], mix[[nm]][i], est[i]))
}
tab <- do.call(rbind, rows)
cat("\n======== summary ========\n")
print(round(tab, 4))
cat(sprintf("\n  mean RMSE across mixtures: %.4f\n", mean(tab[, "RMSE"])))
if (isTRUE(a$keep)) cat("kept workspace: ", tmp, "\n")
