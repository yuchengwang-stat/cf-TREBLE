cf_BLEND_ONE <- function(X, D_X, phi, alpha = 1.001, beta = 1.001,
                    n.iter = 10000, thres = 1e-5){
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  mu <- rep(1/CT, CT)
  psi <- sapply(M, function(x){rep(1/x,x)})
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

cf_BLEND <- function(X, D_X, phi){
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  mu <- rep(1/CT, CT)
  psi <- sapply(M, function(x){rep(1/x,x)})

  weight_SQP <- c(X, (D_X - X))
  L_SQP <- rbind(rlist::list.cbind(phi), rlist::list.cbind(phi_unmethylated))
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
