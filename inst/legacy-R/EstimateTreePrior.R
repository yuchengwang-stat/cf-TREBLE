library(data.table)
append_complements <- function(leaf) {
  n <- length(leaf)
  stopifnot(n >= 1)

  U <- leaf[[n]]
  comps <- lapply(leaf[seq_len(n-1)], function(x) setdiff(U, x))
  leaf <- c(leaf, comps)
  keys <- vapply(leaf, function(v) paste(sort(v), collapse=","), character(1))
  dup_idx <- duplicated(keys)
  for (i in seq_along(leaf)) {
    if (dup_idx[i]) {
      leaf[[i]] <- 0
    }
  }
  return(leaf)
}
#' EstimateTreePrior
#' @description Estimate the prior on the Bayesian tree
#' @export
EstimateTreePrior <- function(BinaryTree, leaf, layer, file_name = "Estimation",tor = 1e-6)
{
  leaf_all <- append_complements(leaf)
  n_node <- length(leaf)
  n_col <- length(leaf_all)
  Loglik_full <- matrix(t(fread(paste0(file_name,"_","ll",".txt"), sep = "\t"))
                             ,ncol = n_col)
  Loglik_truncate <- matrix(t(fread(paste0(file_name,"_","ll.truncate",".txt"), sep = "\t"))
                        ,ncol = n_col)
  n_leaf <- length(leaf[[length(leaf)]])
  n_nonleaf <- n_node - n_leaf
  n_other <- n_col - n_node
  n_1vsrest_tree <- length(leaf_all[!sapply(leaf_all, function(x) identical(x, 0))]) - n_node
  zero_nodes_idx <- which(sapply(leaf_all, function(x) identical(x, 0)))
  pi.t0 <- 0.5
  pi.j0 <- c(rep(1,n_leaf),pi.t0 * rep(1/n_nonleaf,n_nonleaf))
  pi.t <- c(pi.t0,(1-pi.t0) * rep(1/n_1vsrest_tree,n_other))
  pi.t[zero_nodes_idx - n_node + 1] <- 0
  i <- 0
  error <- 999999
  nCpG <- nrow(Loglik_full)
  while(error >= tor)
  {
    i <- i + 1
    p_all <- calculate.p.gj(node.binary, pi.j0, Loglik_full,layer)
    p.t = calculate.p.t(pi.t, pi.j0, n_col, zero_nodes_idx, Loglik_full,Loglik_truncate, p_all[[2]])
    pi.t1 = as.vector(colSums(p.t)/nCpG)
    p = p_all[[1]]
    w <- calculate.w.gj(node.binary, p)
    pi.j1 <- calculate.pi.j(p, w, p.t[,1])
    error <- sum((c(pi.j1,pi.t1) - c(pi.j0,pi.t))^2)
    print(paste0("Number of iteration:",i))
    print(error)
    pi.j0 <- pi.j1
    pi.t <- pi.t1
  }
  saveRDS(c(pi.j0,pi.t),paste0(file_name,"_tree_pi.RDS"))
}

