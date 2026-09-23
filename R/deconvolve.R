## cf-TREBLE pipeline -- deconvolution -----------------------------------------

## .beta is a flat uint8 stream of (methylated, total) pairs, one pair per CpG.
## The CpG count is the file size / 2; pass n_cpg only to assert an expected value.
read_beta <- function(path, n_cpg = NULL) {
  sz <- file.size(path)
  if (is.na(sz)) die("cannot stat ", path)
  if (sz %% 2L != 0L) die(path, " has an odd byte count (", sz, "); not a .beta file")
  n <- as.integer(sz %/% 2L)
  if (!is.null(n_cpg) && !is.na(n_cpg) && n_cpg != n)
    die(path, " holds ", n, " CpGs but the config asserts ", n_cpg)
  raw <- readBin(path, "integer", n = 2L * n, size = 1L, signed = FALSE)
  matrix(raw, ncol = 2L, byrow = TRUE)      # [, 1] = methylated, [, 2] = total
}

## BLEND for methylation.  A cell-free sample is modelled as a mixture over cell
## types, each of which is itself a mixture over the reference subjects of that
## type, so the fitted quantities are a cellular fraction per cell type and a
## mixing proportion over that type's reference panel.
##
##   X     methylated reads of the sample, length G
##   D_X   total reads of the sample, length G
##   phi   one G x M_t beta-value matrix per cell type, named by cell type
##
## `cf_BLEND` poses this as a convex problem over the stacked (methylated,
## unmethylated) references and solves it with mixsqp.  The body is the original
## cf-TREBLE implementation; the only edit is `do.call(cbind, .)` in place of
## `rlist::list.cbind(.)`, which is the same operation without the dependency.
cf_BLEND <- function(X, D_X, phi){
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  # Initialization
  mu <- rep(1/CT, CT)
  psi <- sapply(M, function(x){rep(1/x,x)})

  weight_SQP <- c(X, (D_X - X))
  L_SQP <- rbind(do.call(cbind, phi), do.call(cbind, phi_unmethylated))
  mu_star <- mixsqp::mixsqp(as.matrix(L_SQP), weight_SQP, control = list(verbose = FALSE))
  mu_star <- mu_star$x
  mu_star <- mu_star + 1e-22
  grp <- factor(rep(names(M), times = M),
                levels = names(M))
  mu <- tapply(mu_star, grp, sum)
  mu_list <- split(mu_star, grp)
  psi <- Map(function(x, s) x/s, mu_list, mu)

  names(mu) <- ct.names
  names(psi) <- ct.names
  return(list("cellular_frac"=mu,
              "ref_mixing_prop"=psi))
}

## Build the per-cell-type reference list from a signature table.
## `drop_cols` removes sample columns (e.g. the CVS block) before deconvolving.
signature_to_phi <- function(beta_mat, sample_list, drop_cols = integer(0)) {
  groups <- sample_groups(sample_list)
  keep <- setdiff(seq_along(groups), drop_cols)
  groups <- groups[keep]
  cols <- split(seq_along(groups), groups)
  names(cols) <- names(sample_list)[as.integer(names(cols))]
  lapply(cols, function(ix) beta_mat[, keep[ix], drop = FALSE])
}

## Load a signature for deconvolution.  Two forms are accepted:
##
##   *.rds   what stage `signature` writes: a list with $table (carrying `index`)
##           and $beta (one column per cell type).
##   a table an externally produced signature, as tsv/csv, optionally gzipped:
##           an `index` column giving each CpG's row number in the .beta files
##           to be deconvolved, and one column per cell type holding its beta
##           value in [0, 1].  Metadata columns written by this pipeline are
##           ignored, so signature.txt.gz can be read straight back in.
##
## The table form needs no reference panel and no sample_list -- a signature from
## anywhere can be used, as long as its CpG indexing matches the .beta files.
SIG_META <- c("index", "local_idx", "target", "kind", "n_target", "median_iqr",
              "celltype", "node", "onevsrest", "semi_pair", "score")
read_signature <- function(path, celltypes = NULL) {
  if (!file.exists(path)) die("signature not found: ", path)
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    S <- readRDS(path)
    if (is.null(S$table$index) || is.null(S$beta))
      die(path, " is an .rds but not a signature (needs $table$index and $beta)")
    beta <- as.matrix(S$beta); idx <- as.integer(S$table$index)
  } else {
    d <- data.table::fread(path, showProgress = FALSE, data.table = FALSE)
    if (!"index" %in% names(d))
      die(path, " has no `index` column; it must give each CpG's row number in ",
          "the .beta files (columns found: ", paste(head(names(d), 8), collapse = ", "), ")")
    idx <- as.integer(d$index)
    keep <- setdiff(names(d), SIG_META)
    keep <- keep[!grepl("^(mu|sigma)_", keep)]
    if (!length(keep)) die(path, " has an index column but no cell-type columns")
    beta <- as.matrix(d[, keep, drop = FALSE])
    if (!is.numeric(beta)) die(path, " has non-numeric cell-type columns: ",
                               paste(keep[!vapply(d[keep], is.numeric, TRUE)], collapse = ", "))
  }
  if (!is.null(celltypes)) {
    miss <- setdiff(celltypes, colnames(beta))
    if (length(miss)) die("cell types not in the signature: ", paste(miss, collapse = ", "))
    beta <- beta[, celltypes, drop = FALSE]
  }
  if (anyNA(idx) || any(idx < 1L)) die(path, ": `index` must be positive row numbers")
  bad <- is.finite(beta) & (beta < 0 | beta > 1)
  if (any(bad)) die(path, ": ", sum(bad), " beta values fall outside [0, 1]")
  if (anyNA(beta)) die(path, ": ", sum(is.na(beta)), " missing beta values")
  list(index = idx, beta = beta, celltypes = colnames(beta))
}

## One phi matrix per cell type, one column each.  Within-cell-type reference
## columns only matter for the mixing proportions; replicating a cell type's
## column across its reference samples leaves the cellular fractions unchanged
## (checked to 4e-15 in inst/tests/test_deconvolve.R).
signature_phi <- function(beta) {
  phi <- lapply(seq_len(ncol(beta)), function(k) beta[, k, drop = FALSE])
  setNames(phi, colnames(beta))
}

## Assemble the subject-specific signature written by stage `subject` into one
## marker x reference-sample matrix, aligned to `index`.
##
## BLEND models a cell type as a mixture over its own reference subjects, so it
## wants one column per subject rather than one per cell type; that is what this
## provides.  Stage `subject` runs over the per-chunk marker rows, where a CpG
## selected for two nodes appears twice, but its posterior depends only on the
## CpG, so the duplicates are identical and the first is kept.
##
## A subject with no reads at a CpG has no posterior there and arrives as NA.
## `na` says what to do: "drop" removes those CpGs from the deconvolution (the
## default, and what keeps every column comparable), "celltype" falls back to the
## cell type's own value from the signature for that subject at that CpG.
read_subject <- function(dir, index, beta = NULL, na = c("drop", "celltype")) {
  na <- match.arg(na)
  files <- sort(list.files(dir, "^subject_c[0-9]+[.]rds$", full.names = TRUE))
  if (!length(files)) die("no subject_c*.rds in ", dir,
                          " -- run stage `subject` before deconvolving with --subject")
  parts <- lapply(files, readRDS)
  P   <- do.call(rbind, lapply(parts, `[[`, "p_subj"))
  idx <- unlist(lapply(parts, `[[`, "index"), use.names = FALSE)
  samples <- parts[[1]]$samples
  if (!all(vapply(parts, function(p) identical(p$samples, samples), TRUE)))
    die("the subject chunks in ", dir, " do not agree on the reference samples")
  keep <- !duplicated(idx)
  P <- P[keep, , drop = FALSE]; idx <- idx[keep]
  rows <- match(index, idx)
  if (anyNA(rows))
    die(sum(is.na(rows)), " of ", length(index), " signature CpGs have no subject-specific ",
        "estimate in ", dir, " -- stage `subject` was run on a different marker set")
  P <- P[rows, , drop = FALSE]
  colnames(P) <- samples
  list(beta = P, index = index, samples = samples, na = na)
}

## Split a marker x sample matrix into one matrix per cell type, which is the
## `phi` cf_BLEND takes: entry t is that cell type's reference panel.
## `drop_cols` removes reference samples first, indexed into unlist(sample_list).
subject_phi <- function(P, sample_list, drop_cols = integer(0)) {
  if (ncol(P) != length(unlist(sample_list)))
    die("the subject-specific signature has ", ncol(P), " reference samples but ",
        "sample_list describes ", length(unlist(sample_list)))
  signature_to_phi(P, sample_list, drop_cols)
}

## Deconvolve one sample against `phi`, restricted to the signature CpGs.
deconvolve_beta <- function(beta_path, n_cpg, cpg_index, phi) {
  dat <- read_beta(beta_path, n_cpg)
  if (max(cpg_index) > nrow(dat))
    die(basename(beta_path), " has ", nrow(dat), " CpGs but the signature indexes ",
        max(cpg_index), " -- signature and .beta are on different CpG sets")
  dat <- dat[cpg_index, , drop = FALSE]
  cf_BLEND(dat[, 1], dat[, 2], phi)$cellular_frac
}

## Also carried over from the original cf_BLEND.R, and not used by the pipeline:
## an EM for the same model, kept so that code written against it still runs.
cf_BLEND_em <- function(X, D_X, phi, alpha = 1.001, beta = 1.001,
                        n.iter = 10000, thres = 1e-5) {
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  # Initialization
  mu <- rep(1/CT, CT)
  # lapply, not sapply: when every cell type has the same number of reference
  # subjects sapply simplifies to a matrix, and psi[[t]] is then one scalar
  # rather than that type's vector.  The panels this was written against had
  # unequal sizes, so the simplification never fired there.
  psi <- lapply(M, function(x){rep(1/x,x)})
  U_1 <- list()
  U_0 <- list()
  for(t in 1:CT){
    U_1 <- c(U_1, list(matrix(0, nrow = G, ncol = M[t])))
    U_0 <- c(U_0, list(matrix(0, nrow = G, ncol = M[t])))
  }
  n.converge <- 0 # record number of iterations till convergence

  for(i in 1:n.iter){
    if((i %% 200)==0){gc()}
    n.converge <- n.converge + 1
    # Update U_1 and U_0
    for(t in 1:CT){
      U_1[[t]] <- mu[t]*t(apply(phi[[t]], 1, function(x){x*psi[[t]]}))
      U_0[[t]] <- mu[t]*t(apply(phi_unmethylated[[t]], 1, function(x){x*psi[[t]]}))
    }
    U_1_sum <- Reduce("+",lapply(U_1, rowSums))
    tmp_1 <- X/U_1_sum
    U_0_sum <- Reduce("+",lapply(U_0, rowSums))
    tmp_0 <- (D_X - X)/U_0_sum
    for(t in 1:CT){
      U_1[[t]] <- apply(U_1[[t]], 2, function(x){x*tmp_1})
      U_0[[t]] <- apply(U_0[[t]], 2, function(x){x*tmp_0})
    }
    U <- U_1
    for(t in 1:CT){
      U[[t]] <- U[[t]] + U_0[[t]]
    }
    # Record the last mu before updating it
    mu_old <- mu
    # Update mu
    mu <- unlist(lapply(U, sum)) + (alpha - 1)
    mu <- mu/sum(mu)
    # Update psi
    for(t in 1:CT){
      psi[[t]] <- colSums(U[[t]]) + (beta - 1)
      psi[[t]] <- psi[[t]]/sum(psi[[t]])
    }
    if(sum(abs(mu - mu_old))<thres){
      names(mu) <- ct.names
      names(psi) <- ct.names
      return(list("cellular_frac"=mu,
                  "ref_mixing_prop"=psi,
                  "iterations"=n.converge))
    }
  }
  names(mu) <- ct.names
  names(psi) <- ct.names
  gc()
  return(list("cellular_frac"=mu,
              "ref_mixing_prop"=psi,
              "iterations"=n.converge))
}

## The name this solver had in the original scripts, kept so existing code runs.
cf_BLEND_ONE <- cf_BLEND_em
