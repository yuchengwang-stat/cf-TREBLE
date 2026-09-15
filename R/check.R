## cfTREBLE -- validate a configuration's inputs before spending any compute ----
##
## The expensive stage runs for a day per chunk.  Every mistake this catches is
## one that would otherwise surface late, or not at all: a tree whose leaves are
## not the first n_leaf nodes still "works" and produces a subtly wrong
## signature.  So the checks are strict and say what is wrong, not just that
## something is.

.ck <- function(state, ok, what, detail = "") {
  state$rows[[length(state$rows) + 1L]] <-
    data.frame(check = what, result = if (isTRUE(ok)) "ok" else "FAIL",
               detail = detail, stringsAsFactors = FALSE)
  if (!isTRUE(ok)) state$bad <- state$bad + 1L
  invisible(state)
}

#' Validate the inputs a configuration points at
#'
#' @param cfg a configuration list from `load_config`
#' @param chunks which count-matrix chunks to open (default: the first)
#' @param deep also read every configured chunk, not just `chunks`
#' @return invisibly, a data frame of checks; prints a report
#' @export
check_inputs <- function(cfg, chunks = NULL, deep = FALSE) {
  st <- new.env(); st$rows <- list(); st$bad <- 0L
  ex <- function(p) !is.null(p) && nzchar(p) && file.exists(p)

  ## ---- sample list ------------------------------------------------------------
  slf <- cfg$data$sample_list
  .ck(st, ex(slf), "sample_list file exists", slf)
  sl <- NULL
  if (ex(slf)) {
    sl <- readRDS(slf)
    .ck(st, is.list(sl) && !is.null(names(sl)) && all(nzchar(names(sl))),
        "sample_list is a named list", sprintf("%d entries", length(sl)))
    .ck(st, all(lengths(sl) >= 1), "every cell type has >= 1 sample",
        sprintf("min %d, max %d", min(lengths(sl)), max(lengths(sl))))
    .ck(st, !anyDuplicated(names(sl)), "cell type names are unique")
    .ck(st, !anyDuplicated(unlist(sl)), "sample names are unique",
        sprintf("%d samples total", length(unlist(sl))))
  }

  ## ---- tree -------------------------------------------------------------------
  fl <- cfg$tree$from_lists
  tr <- NULL
  if (ex(fl$children)) {
    ch <- readRDS(fl$children)
    pa <- if (ex(fl$parent)) readRDS(fl$parent) else NULL
    n <- length(ch)
    leaves <- which(lengths(ch) == 1L)
    .ck(st, identical(sort(leaves), seq_len(length(leaves))),
        "leaves are nodes 1..n_leaf",
        if (identical(sort(leaves), seq_len(length(leaves)))) ""
        else sprintf("leaf nodes are %s", paste(head(leaves, 8), collapse = ",")))
    .ck(st, n == 2L * length(leaves) - 1L, "tree is binary",
        sprintf("%d nodes, %d leaves; expected %d nodes", n, length(leaves),
                2L * length(leaves) - 1L))
    if (!is.null(sl))
      .ck(st, length(leaves) == length(sl),
          "one leaf per cell type",
          sprintf("%d leaves, %d cell types", length(leaves), length(sl)))
    if (!is.null(pa))
      .ck(st, as.integer(unlist(pa))[n] == n, "root is its own parent",
          sprintf("parent[%d] = %s", n, as.integer(unlist(pa))[n]))
    tr <- tryCatch(tree_from_lists(fl$children, fl$parent, fl$celltypes),
                   error = function(e) conditionMessage(e))
    .ck(st, !is.character(tr), "tree builds and validates",
        if (is.character(tr)) tr else
          sprintf("%d leaves, %d nodes, %d likelihood columns, %d layers",
                  tr$n_leaf, tr$n_node, tr$n_col, length(tr$layer)))
    if (is.character(tr)) tr <- NULL
  } else {
    .ck(st, FALSE, "tree children file exists", fl$children %||% "(unset)")
  }

  ## ---- priors -----------------------------------------------------------------
  if (ex(cfg$priors$mu)) {
    mp <- readRDS(cfg$priors$mu)
    have <- all(c("mu", "sigma", "lambda") %in% names(mp))
    .ck(st, have, "mu prior has $mu, $sigma, $lambda",
        paste(names(mp), collapse = ", "))
    if (have) {
      k <- length(mp$mu)
      .ck(st, length(mp$sigma) == k && length(mp$lambda) == k,
          "mu prior components agree in length", sprintf("%d components", k))
      .ck(st, abs(sum(mp$lambda) - 1) < 1e-6, "mu prior weights sum to 1",
          sprintf("sum = %.6f", sum(mp$lambda)))
      .ck(st, all(mp$sigma > 0), "mu prior sigmas are positive")
    }
  } else .ck(st, FALSE, "mu prior file exists", cfg$priors$mu %||% "(unset)")

  if (ex(cfg$priors$sigma)) {
    sp <- readRDS(cfg$priors$sigma)
    have <- all(c("grid", "pi") %in% names(sp))
    .ck(st, have, "sigma prior has $grid and $pi", paste(names(sp), collapse = ", "))
    if (have) {
      .ck(st, length(sp$grid) == length(sp$pi), "grid and pi agree in length",
          sprintf("%d points", length(sp$grid)))
      .ck(st, !is.unsorted(sp$grid), "grid is ascending",
          paste(signif(sp$grid, 3), collapse = " "))
      .ck(st, sp$grid[1] == 0, "grid starts at 0",
          "the zero-variance component is what the truncated likelihood uses")
      .ck(st, abs(sum(sp$pi) - 1) < 1e-6, "sigma prior weights sum to 1",
          sprintf("sum = %.6f", sum(sp$pi)))
      tr_i <- cfg$priors$sigma_truncate
      if (is.null(tr_i) || identical(tr_i, "auto"))
        .ck(st, any(sp$grid < 0.5), "auto truncation selects something",
            sprintf("%d of %d grid points below 0.5", sum(sp$grid < 0.5), length(sp$grid)))
    }
  } else .ck(st, FALSE, "sigma prior file exists", cfg$priors$sigma %||% "(unset)")

  ## ---- count matrices ----------------------------------------------------------
  all_chunks <- as.integer(unlist(cfg$data$chunks))
  if (is.null(chunks)) chunks <- if (deep) all_chunks else all_chunks[1]
  n_samp <- if (!is.null(sl)) length(unlist(sl)) else NA_integer_
  heights <- integer(0)
  for (k in chunks) {
    mf <- sprintf(cfg$data$m_chunks, k); nf <- sprintf(cfg$data$n_chunks, k)
    if (!ex(mf) || !ex(nf)) { .ck(st, FALSE, sprintf("chunk %d present", k),
                                  paste(basename(c(mf, nf)), collapse = " / ")); next }
    M <- readRDS(mf); N <- readRDS(nf)
    heights <- c(heights, nrow(M))
    .ck(st, identical(dim(M), dim(N)), sprintf("chunk %d: M and N same shape", k),
        sprintf("%d x %d", nrow(M), ncol(M)))
    if (!is.na(n_samp))
      .ck(st, ncol(M) == n_samp, sprintf("chunk %d: columns match sample_list", k),
          sprintf("%d columns, sample_list implies %d", ncol(M), n_samp))
    .ck(st, all(M <= N, na.rm = TRUE), sprintf("chunk %d: M <= N everywhere", k))
    .ck(st, !anyNA(M) && !anyNA(N), sprintf("chunk %d: no NA", k))
    .ck(st, min(M, na.rm = TRUE) >= 0, sprintf("chunk %d: counts non-negative", k))
    rm(M, N); gc()
  }
  stride <- as.integer(cfg$data$chunk_stride %||% 0L)
  if (length(heights) && stride)
    .ck(st, all(heights <= stride), "chunk heights fit chunk_stride",
        sprintf("heights %s, stride %d", paste(unique(heights), collapse = ","), stride))

  ## ---- report ------------------------------------------------------------------
  out <- do.call(rbind, st$rows)
  cat("\n")
  for (i in seq_len(nrow(out)))
    cat(sprintf("  %-4s %-42s %s\n",
                if (out$result[i] == "ok") "ok" else "FAIL",
                out$check[i], out$detail[i]))
  cat(sprintf("\n  %d checks, %d failed\n", nrow(out), st$bad))
  if (st$bad) cat("  Fix these before running anything expensive.\n")
  invisible(out)
}
