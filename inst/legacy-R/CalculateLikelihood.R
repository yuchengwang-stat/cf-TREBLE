suppressPackageStartupMessages({
  library(Matrix); library(expm); library(matrixcalc); library(numDeriv)
  library(pbv);    library(fastGHQuad); library(pracma);
  library(data.table); library(parallel); library(RhpcBLASctl)
})

#' CalculateLikelihood
#' @description Estimate the prior on the Bayesian tree
#' @export
CalculateLikelihood <- function(M, N, sample_list, leaf, mu.prior, sigma.prior, ncores = 1,
                                GH.point = 5, file_name = "Estimation")
{
  Leaf.tree <- leaf
  samples <- sample_list
  samples_num <- lengths(samples)
  names(samples_num) <- NULL
  samples_num <- unlist(samples_num)
  samples_num <- rep(1:length(samples_num), times = samples_num)
  split_cols  <- split(seq_along(samples_num), samples_num)
  delta       <- mu.prior$mu
  tau2        <- mu.prior$sigma^2
  lambda      <- mu.prior$lambda
  sigma.grid  <- sigma.prior$grid
  sigma.pi    <- sigma.prior$pi
  Gauss.grid  <- gaussHermite(5)$x
  Gauss.w     <- gaussHermite(5)$w
  fields <- c("ll","ll.truncate","mu","mu2","sigma",
                   "p","p2","p0","pg0","alpha","pi.Inf.alpha")
  safe_Post <- function(Mrow, Urow) {
    tryCatch(
      PostExp_fast(M=Mrow, U=Urow,
                   sigma.grid=sigma.grid, sigma.pi=sigma.pi,
                   mu.mu=delta, mu.tau2=tau2, mu.pi=lambda,
                   Gauss.grid=Gauss.grid, Gauss.w=Gauss.w,
                   ind.sigma.truncate=1:6),
      error=function(e)
        list(ll=rep(NA_real_,1), ll.truncate=rep(NA_real_,1),
             mu=NA, mu2=NA, sigma=NA, p=NA, p2=NA,
             p0=NA, pg0=NA, alpha=NA, pi.Inf.alpha=NA)
    )
  }
  append_complements <- function(leaf) {
    n <- length(leaf)
    stopifnot(n >= 1)

    U <- leaf[[n]]
    comps <- lapply(leaf[seq_len(n-1)], function(x) setdiff(U, x))
    leaf <- c(leaf, comps)

    return(leaf)
  }
  leaf <- append_complements(leaf)
  n <- nrow(M)
  chunk_size  <- ceiling(n/ncores)
  flush_n     <- 600000


  buf   <- vector("list", length(fields)); names(buf) <- fields
  add2buf <- function(res) {
    for (nm in fields) buf[[nm]][[length(buf[[nm]])+1]] <<- res[[nm]]
  }
  flush_buf <- function() {
    if (length(buf[[1]])==0) return()
    for (nm in fields) {
      dt <- transpose(buf[[nm]])
      fwrite(dt, file=sprintf("%s_%s.txt", file_name, nm),
             append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
      buf[[nm]] <<- list()
    }
    print("flushed!")
  }
  n_in_buf <- 0
  U <- N - M
  M_ls  <- lapply(split_cols, \(ix) M[, ix, drop=FALSE])
  U_ls  <- lapply(split_cols, \(ix) U[, ix, drop=FALSE])
  for (i in seq_along(Leaf.tree)) {
    leaf_timer <- Sys.time()
    idx   <- Leaf.tree[[i]]
    M_now <- do.call(cbind, M_ls[idx])
    U_now <- do.call(cbind, U_ls[idx])

    rows  <- split(seq_len(nrow(M_now)),
                   ceiling(seq_len(nrow(M_now))/chunk_size))

    res   <- mclapply(
      rows,
      \(row_idx) {
        apply(cbind(M_now[row_idx, , drop=FALSE],
                    U_now[row_idx, , drop=FALSE]), 1,
              \(ru) safe_Post(ru[seq_len(ncol(M_now))],
                              ru[-seq_len(ncol(M_now))]))
      },
      mc.cores=ncores, mc.preschedule=FALSE
    )
    chunk_lengths <- sapply(res, length)
    for (chunk in res) {
      lapply(chunk, add2buf)
      n_in_buf <- n_in_buf + length(chunk)
      if (n_in_buf >= flush_n) { flush_buf(); n_in_buf <- 0 }
    }
    if (n_in_buf >= flush_n) { flush_buf(); n_in_buf <- 0 }
    cat("leaf", i,
        " rows=", nrow(M_now),
        " time=", round(difftime(Sys.time(), leaf_timer, units="secs"),1), "s\n")
  }

  flush_buf()
}

#' CalculateIndp
#' @description Estimate the prior on the Bayesian tree
#' @export
CalculateIndp <- function(M, N, index, sample_list, leaf, mu.prior, sigma.prior, ncores = 1,
                                GH.point = 5, file_name = "Estimation")
{
  Leaf.tree <- leaf
  samples <- sample_list
  samples_num <- lengths(samples)
  names(samples_num) <- NULL
  samples_num <- unlist(samples_num)
  samples_num <- rep(1:length(samples_num), times = samples_num)
  split_cols  <- split(seq_along(samples_num), samples_num)
  delta       <- mu.prior$mu
  tau2        <- mu.prior$sigma^2
  lambda      <- mu.prior$lambda
  sigma.grid  <- sigma.prior$grid
  sigma.pi    <- sigma.prior$pi
  Gauss.grid  <- gaussHermite(5)$x
  Gauss.w     <- gaussHermite(5)$w
    fields <- c("Ind")
    safe_Post <- function(Mrow, Urow) {
      tryCatch(
        PostExp_IndSig(
          M = Mrow,
          U = Urow,
          sigma.grid = sigma.grid,
          sigma.pi   = sigma.pi,
          mu.mu      = delta,
          mu.tau2    = tau2,
          mu.pi      = lambda,
          Gauss.grid = Gauss.grid,
          Gauss.w    = Gauss.w
        ),
        error = function(e) rep(NA_real_, length(Mrow))
      )
    }
  append_complements <- function(leaf) {
    n <- length(leaf)
    stopifnot(n >= 1)

    U <- leaf[[n]]
    comps <- lapply(leaf[seq_len(n-1)], function(x) setdiff(U, x))
    leaf <- c(leaf, comps)

    return(leaf)
  }
  leaf <- append_complements(leaf)
  n <- nrow(M)
  chunk_size  <- ceiling(n/ncores)
  flush_n     <- 600000
  buf   <- vector("list", length(fields)); names(buf) <- fields
  add2buf <- function(res) {
    for (nm in fields) buf[[nm]][[length(buf[[nm]])+1]] <<- res[[nm]]
  }
  flush_buf <- function() {
    if (length(buf[[1]])==0) return()
    for (nm in fields) {
      dt <- transpose(buf[[nm]])
      fwrite(dt, file=sprintf("%s_%s.txt", file_name, nm),
             append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
      buf[[nm]] <<- list()
    }
    print("flushed!")
  }
  n_in_buf <- 0
  M <- M[index,]
  N <- N[index,]
  U <- N - M
  M_ls  <- lapply(split_cols, \(ix) M[, ix, drop=FALSE])
  U_ls  <- lapply(split_cols, \(ix) U[, ix, drop=FALSE])
  for (i in seq_along(Leaf.tree)) {
    leaf_timer <- Sys.time()
    idx   <- Leaf.tree[[i]]
    M_now <- do.call(cbind, M_ls[idx])
    U_now <- do.call(cbind, U_ls[idx])

    rows  <- split(seq_len(nrow(M_now)),
                   ceiling(seq_len(nrow(M_now))/chunk_size))

    res   <- mclapply(
      rows,
      \(row_idx) {
        apply(cbind(M_now[row_idx, , drop=FALSE],
                    U_now[row_idx, , drop=FALSE]), 1,
              \(ru) safe_Post(ru[seq_len(ncol(M_now))],
                              ru[-seq_len(ncol(M_now))]))
      },
      mc.cores=ncores, mc.preschedule=FALSE
    )
    chunk_lengths <- sapply(res, length)
    for (chunk in res) {
      lapply(chunk, add2buf)
      n_in_buf <- n_in_buf + length(chunk)
      if (n_in_buf >= flush_n) { flush_buf(); n_in_buf <- 0 }
    }
    if (n_in_buf >= flush_n) { flush_buf(); n_in_buf <- 0 }
    cat("leaf", i,
        " rows=", nrow(M_now),
        " time=", round(difftime(Sys.time(), leaf_timer, units="secs"),1), "s\n")
  }

  flush_buf()
}
