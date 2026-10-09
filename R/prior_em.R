calculate.p.gj <- function(children, pi.j, log.lik, layer) {
  split_ll <- log.lik
  inductive <- function(index, lls) {
    kid <- children[[index]]
    if (length(kid) == 1L) return(lls[, index])
    l <- kid[1]; r <- kid[2]
    summand <- cbind(lls[, l]     + lls[, r]     + log(1 - pi.j[l]) + log(1 - pi.j[r]),
                     lls[, l]     + log.lik[, r] + log(1 - pi.j[l]) + log(pi.j[r]),
                     log.lik[, l] + lls[, r]     + log(pi.j[l])     + log(1 - pi.j[r]),
                     log.lik[, l] + log.lik[, r] + log(pi.j[l])     + log(pi.j[r]))
    C <- apply(summand, 1, max)
    log(rowSums(exp(summand - C))) + C
  }
  for (lay in layer) split_ll[, lay] <- mapply(inductive, lay, MoreArgs = list(split_ll))
  n_node <- length(children)
  out <- sweep((split_ll - log.lik)[, seq_len(n_node)], 2, log(1 - pi.j) - log(pi.j), "+")
  list(1 / (1 + exp(out)), split_ll[, n_node])
}

calculate.p.t <- function(pi.t, pi.z, tree, ll_full, ll_trunc, ll_root_1) {
  n_node <- tree$n_node; k <- tree$comp_of
  ll_Tk <- ll_trunc[, n_node + seq_along(k), drop = FALSE] + ll_trunc[, k, drop = FALSE]
  t1 <- ll_full[, n_node] + log(pi.z[n_node])
  t2 <- ll_root_1        + log(1 - pi.z[n_node])
  C  <- pmax(t1, t2)
  ll_T0 <- C + log(exp(t1 - C) + exp(t2 - C))
  un <- cbind(ll_T0 + log(pi.t[1]), sweep(ll_Tk, 2, log(pi.t[-1]), "+"))
  CC <- apply(un, 1, max)
  exp(un - (log(rowSums(exp(un - CC))) + CC))
}

calculate.w.gj <- function(children, p) {
  out <- 1 - p
  recurse <- function(index, acc) {
    if (length(children[[index]]) == 1L) { out[, index] <<- acc; return(invisible()) }
    nxt <- acc * out[, index]
    out[, index] <<- acc
    for (k in children[[index]]) recurse(k, nxt)
  }
  recurse(length(children), rep(1, nrow(out)))
  out
}

calculate.pi.j <- function(p, w, p.t0) colSums(p.t0 * p * w) / colSums(p.t0 * w)

estimate_tree_prior <- function(tree, ll_full, ll_trunc, tol = 1e-6, max_iter = 500) {
  n_node <- tree$n_node; n_leaf <- tree$n_leaf
  n_internal <- n_node - n_leaf
  n_comp <- length(tree$comp_of)
  if (ncol(ll_full) != tree$n_col)
    die("likelihood has ", ncol(ll_full), " columns but the tree needs ", tree$n_col)

  pi.t0 <- 0.5
  pi.j <- c(rep(1, n_leaf), pi.t0 * rep(1 / n_internal, n_internal))
  pi.t <- c(pi.t0, (1 - pi.t0) * rep(1 / n_comp, n_comp))
  nCpG <- nrow(ll_full)
  trace <- numeric(0)
  eps <- 1e-10
  internal <- seq.int(n_leaf + 1L, n_node)
  clamp <- function(v) pmin(pmax(v, eps), 1 - eps)

  for (it in seq_len(max_iter)) {
    pa   <- calculate.p.gj(tree$children, pi.j, ll_full, tree$layer)
    p.t  <- calculate.p.t(pi.t, pi.j, tree, ll_full, ll_trunc, pa[[2]])
    pi.t.new <- as.vector(colSums(p.t) / nCpG)
    w    <- calculate.w.gj(tree$children, pa[[1]])
    pi.j.new <- calculate.pi.j(pa[[1]], w, p.t[, 1])
    pi.j.new[internal] <- clamp(pi.j.new[internal])
    pi.t.new <- pmax(pi.t.new, eps); pi.t.new <- pi.t.new / sum(pi.t.new)
    err <- sum((c(pi.j.new, pi.t.new) - c(pi.j, pi.t))^2)
    if (!is.finite(err)) {
      bad <- c(which(!is.finite(pi.j.new)), n_node + which(!is.finite(pi.t.new)))
      die("EM diverged at iteration ", it,
          " -- non-finite pi at position(s) ", paste(head(bad, 10), collapse = ","),
          ". Usually too few CpGs: raise --max-rows.")
    }
    trace <- c(trace, err)
    log_msg(sprintf("  EM iter %3d   err = %.3e", it, err))
    pi.j <- pi.j.new; pi.t <- pi.t.new
    if (err < tol) break
  }
  if (tail(trace, 1) >= tol)
    log_msg("WARNING: EM hit max_iter (", max_iter, ") without reaching tol")
  list(pi.j = pi.j, pi.t = pi.t, trace = trace, converged = tail(trace, 1) < tol,
       n_cpg = nCpG, n_node = n_node)
}

estimate_tree_prior_streaming <- function(tree, chunks, loader, tol = 1e-6,
                                          max_iter = 500, cache = FALSE,
                                          progress = TRUE) {
  n_node <- tree$n_node; n_leaf <- tree$n_leaf
  n_internal <- n_node - n_leaf; n_comp <- length(tree$comp_of)
  eps <- 1e-10; internal <- seq.int(n_leaf + 1L, n_node)
  clamp <- function(v) pmin(pmax(v, eps), 1 - eps)

  pi.t0 <- 0.5
  pi.j <- c(rep(1, n_leaf), pi.t0 * rep(1 / n_internal, n_internal))
  pi.t <- c(pi.t0, (1 - pi.t0) * rep(1 / n_comp, n_comp))
  store <- if (cache) vector("list", length(chunks)) else NULL
  trace <- numeric(0); nCpG_total <- NA_integer_

  for (it in seq_len(max_iter)) {
    acc_t <- numeric(1L + n_comp); acc_num <- numeric(n_node)
    acc_den <- numeric(n_node); nC <- 0L
    t0 <- Sys.time()
    for (ci in seq_along(chunks)) {
      k <- chunks[ci]
      dat <- if (cache && !is.null(store[[ci]])) store[[ci]] else loader(k)
      if (is.null(dat)) next
      if (cache && is.null(store[[ci]])) store[[ci]] <- dat
      pa  <- calculate.p.gj(tree$children, pi.j, dat$ll, tree$layer)
      p.t <- calculate.p.t(pi.t, pi.j, tree, dat$ll, dat$llt, pa[[2]])
      w   <- calculate.w.gj(tree$children, pa[[1]])
      acc_t   <- acc_t   + colSums(p.t)
      acc_num <- acc_num + colSums(p.t[, 1] * pa[[1]] * w)
      acc_den <- acc_den + colSums(p.t[, 1] * w)
      nC <- nC + nrow(dat$ll)
      rm(pa, p.t, w); if (!cache) rm(dat)
    }
    if (!nC) die("streaming EM saw no data")
    nCpG_total <- nC
    pi.t.new <- acc_t / nC
    pi.j.new <- acc_num / acc_den
    pi.j.new[internal] <- clamp(pi.j.new[internal])
    pi.t.new <- pmax(pi.t.new, eps); pi.t.new <- pi.t.new / sum(pi.t.new)
    err <- sum((c(pi.j.new, pi.t.new) - c(pi.j, pi.t))^2)
    if (!is.finite(err)) die("streaming EM diverged at iteration ", it)
    trace <- c(trace, err)
    if (progress)
      log_msg(sprintf("  EM iter %3d  err = %.3e  (%d CpGs, %.1f s)", it, err, nC,
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    pi.j <- pi.j.new; pi.t <- pi.t.new
    if (err < tol) break
  }
  if (tail(trace, 1) >= tol)
    log_msg("WARNING: EM hit max_iter (", max_iter, ") without reaching tol")
  list(pi.j = pi.j, pi.t = pi.t, trace = trace, converged = tail(trace, 1) < tol,
       n_cpg = nCpG_total, n_node = n_node, streaming = TRUE)
}
