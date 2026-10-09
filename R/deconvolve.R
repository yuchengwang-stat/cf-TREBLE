read_beta <- function(path, n_cpg = NULL) {
  sz <- file.size(path)
  if (is.na(sz)) die("cannot stat ", path)
  if (sz %% 2L != 0L) die(path, " has an odd byte count (", sz, "); not a .beta file")
  n <- as.integer(sz %/% 2L)
  if (!is.null(n_cpg) && !is.na(n_cpg) && n_cpg != n)
    die(path, " holds ", n, " CpGs but the config asserts ", n_cpg)
  raw <- readBin(path, "integer", n = 2L * n, size = 1L, signed = FALSE)
  matrix(raw, ncol = 2L, byrow = TRUE)
}

cf_BLEND <- function(X, D_X, phi){
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

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

signature_to_phi <- function(beta_mat, sample_list, drop_cols = integer(0)) {
  groups <- sample_groups(sample_list)
  keep <- setdiff(seq_along(groups), drop_cols)
  groups <- groups[keep]
  cols <- split(seq_along(groups), groups)
  names(cols) <- names(sample_list)[as.integer(names(cols))]
  lapply(cols, function(ix) beta_mat[, keep[ix], drop = FALSE])
}

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

signature_phi <- function(beta) {
  phi <- lapply(seq_len(ncol(beta)), function(k) beta[, k, drop = FALSE])
  setNames(phi, colnames(beta))
}

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

subject_phi <- function(P, sample_list, drop_cols = integer(0)) {
  if (ncol(P) != length(unlist(sample_list)))
    die("the subject-specific signature has ", ncol(P), " reference samples but ",
        "sample_list describes ", length(unlist(sample_list)))
  signature_to_phi(P, sample_list, drop_cols)
}

deconvolve_beta <- function(beta_path, n_cpg, cpg_index, phi) {
  dat <- read_beta(beta_path, n_cpg)
  if (max(cpg_index) > nrow(dat))
    die(basename(beta_path), " has ", nrow(dat), " CpGs but the signature indexes ",
        max(cpg_index), " -- signature and .beta are on different CpG sets")
  dat <- dat[cpg_index, , drop = FALSE]
  cf_BLEND(dat[, 1], dat[, 2], phi)$cellular_frac
}

cf_BLEND_em <- function(X, D_X, phi, alpha = 1.001, beta = 1.001,
                        n.iter = 10000, thres = 1e-5) {
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  mu <- rep(1/CT, CT)
  psi <- lapply(M, function(x){rep(1/x,x)})
  U_1 <- list()
  U_0 <- list()
  for(t in 1:CT){
    U_1 <- c(U_1, list(matrix(0, nrow = G, ncol = M[t])))
    U_0 <- c(U_0, list(matrix(0, nrow = G, ncol = M[t])))
  }
  n.converge <- 0

  for(i in 1:n.iter){
    if((i %% 200)==0){gc()}
    n.converge <- n.converge + 1
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
    mu_old <- mu
    mu <- unlist(lapply(U, sum)) + (alpha - 1)
    mu <- mu/sum(mu)
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

cf_BLEND_ONE <- cf_BLEND_em
