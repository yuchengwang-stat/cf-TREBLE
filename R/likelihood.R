## cf-TREBLE pipeline -- per-node log-likelihood over CpG chunks ---------------
##
## Output is binary and carries its CpG index, so row order is data rather than a
## convention shared between stages.  One file per (chunk, batch), so a rerun
## skips what is already finished.  The field set is configurable and defaults to
## the four fields the later stages read.

LIK_FIELDS_ALL <- c("ll","ll.truncate","mu","mu2","sigma",
                    "p","p2","p0","pg0","alpha","pi.Inf.alpha")
LIK_FIELDS_MIN <- c("ll","ll.truncate","mu","sigma")

## priors -> the argument bundle PostExp_fast wants
load_priors <- function(cfg) {
  mp <- readRDS(cfg$priors$mu)
  sp <- readRDS(cfg$priors$sigma)
  trunc <- cfg$priors$sigma_truncate
  ind <- if (is.null(trunc) || identical(trunc, "auto")) which(sp$grid < 0.5) else as.integer(unlist(trunc))
  gh <- pracma::gaussHermite(as.integer(cfg$priors$gh_points %||% 5))
  list(delta = mp$mu, tau2 = mp$sigma^2, lambda = mp$lambda,
       sigma.grid = sp$grid, sigma.pi = sp$pi,
       ind.sigma.truncate = ind, gh.x = gh$x, gh.w = gh$w)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## one CpG: posterior summaries under every node of the tree
lik_one <- function(Mrow, Urow, pr, fields) {
  na <- setNames(rep(NA_real_, length(fields)), fields)
  res <- tryCatch(
    PostExp_fast(M = Mrow, U = Urow,
                 sigma.grid = pr$sigma.grid, sigma.pi = pr$sigma.pi,
                 mu.mu = pr$delta, mu.tau2 = pr$tau2, mu.pi = pr$lambda,
                 Gauss.grid = pr$gh.x, Gauss.w = pr$gh.w,
                 ind.sigma.truncate = pr$ind.sigma.truncate),
    error = function(e) NULL)
  if (is.null(res)) return(na)
  vapply(fields, function(f) as.numeric(res[[f]])[1], 0.0)
}

## One (chunk, batch): rows x nodes matrices, one per field.
lik_batch <- function(M, N, tree, groups, pr, fields, ncores, block_size) {
  U <- N - M
  gcols <- split(seq_along(groups), groups)
  M_ls <- lapply(gcols, function(ix) M[, ix, drop = FALSE])
  U_ls <- lapply(gcols, function(ix) U[, ix, drop = FALSE])
  nrows <- nrow(M)
  ## Never leave cores idle: with the configured block size a short batch can
  ## split into fewer blocks than there are cores.
  bs_eff <- max(1L, min(as.integer(block_size), ceiling(nrows / ncores)))
  if (bs_eff != block_size)
    log_msg("block_size ", block_size, " -> ", bs_eff, " so ", ncores,
            " cores each get work (", ceiling(nrows / bs_eff), " blocks)")
  out <- lapply(fields, function(...) matrix(NA_real_, nrows, tree$n_col))
  names(out) <- fields
  t_start <- Sys.time()

  for (j in seq_len(tree$n_col)) {
    idx   <- tree$leaf_all[[j]]
    M_now <- do.call(cbind, M_ls[idx])
    U_now <- do.call(cbind, U_ls[idx])
    blocks <- split(seq_len(nrows), ceiling(seq_len(nrows) / bs_eff))
    res <- parallel::mclapply(blocks, function(rows) {
      vapply(rows, function(r) lik_one(M_now[r, ], U_now[r, ], pr, fields),
             numeric(length(fields)))
    }, mc.cores = ncores, mc.preschedule = FALSE)
    bad <- vapply(res, function(x) inherits(x, "try-error") || !is.matrix(x), TRUE)
    if (any(bad)) die("node ", j, ": ", sum(bad), " block(s) failed in mclapply")
    got <- do.call(cbind, res)                       # fields x nrows
    for (f in seq_along(fields)) out[[f]][, j] <- got[f, ]
    if (j %% 20L == 0L || j == tree$n_col) {
      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      log_msg(sprintf("    node %d/%d  %.1f min elapsed, ~%.1f min left",
                      j, tree$n_col, el, el / j * (tree$n_col - j)))
    }
  }
  out
}

## Filenames are derived, never typed twice.
lik_path <- function(dir, chunk, batch) file.path(dir, sprintf("lik_c%03d_b%03d.rds", chunk, batch))

## Global CpG ids: chunks are laid end to end in the order given by
## config$data$chunks.  A chunk records its own full height when stage
## `likelihood` runs; for chunks that have not run yet we fall back to the
## configured stride, so a single-chunk run does not need all the others.
chunk_offsets <- function(dir, chunks, stride = 900000L) {
  n <- vapply(chunks, function(k) {
    f <- file.path(dir, sprintf("rows_c%03d.txt", k))
    if (file.exists(f)) as.integer(readLines(f)[1]) else NA_integer_
  }, 0L)
  if (anyNA(n)) {
    log_msg("chunk heights unknown for chunk(s) ", paste(chunks[is.na(n)], collapse = ","),
            "; assuming stride ", stride)
    n[is.na(n)] <- as.integer(stride)
  }
  setNames(cumsum(c(0L, head(n, -1))), as.character(chunks))
}

## Load every finished batch of one chunk, row-bound, with global CpG ids attached.
load_chunk <- function(dir, chunk, offset = 0L, fields = NULL, max_rows = 0L) {
  fs <- sort(list.files(dir, pattern = sprintf("^lik_c%03d_b[0-9]+\\.rds$", chunk), full.names = TRUE))
  if (!length(fs)) die("no likelihood output for chunk ", chunk, " in ", dir)
  parts <- lapply(fs, readRDS)
  if (is.null(fields)) fields <- names(parts[[1]]$fields)
  out <- lapply(fields, function(f) do.call(rbind, lapply(parts, function(p) p$fields[[f]])))
  names(out) <- fields
  idx <- offset + unlist(lapply(parts, function(p) p$local_rows), use.names = FALSE)
  if (max_rows > 0L && max_rows < length(idx)) {
    keep <- seq_len(max_rows)
    out <- lapply(out, function(m) m[keep, , drop = FALSE]); idx <- idx[keep]
  }
  list(fields = out, cpg_index = idx, local_rows = idx - offset)
}
