# BLEND-M
# Served as the DNA methylation version of RNA deconvolution method BLEND
# Admixture sample: WGBS DNAm data
# Reference panels: array-based DNAm measurements (beta-values)

## Data Input
### X: mathylated reads of a subject, length G
### D_X: total reads of a subject, length G
### alpha: cellular fraction prior
### beta: reference mixing proportion prior
### phi: list of cell type specific DNAm references, length CT.
###      names of this list represent cell types.
###      Each element of the list is a numerical beta value matrix G X M_t.
### n.iter: maximum number of iterations allowed.
### thres: threshold for convergence identification
###

## Algorithm parameters
### mu: cellular fraction vector, length CT
### psi: list of reference mixing proportions, length CT. each vector has length M_t
### U_1: list, mathylated latent variable, length CT, each element G X M_t
### U_0: list, unmathylated latent variable, length CT, each element G X M_t


cf_BLEND_ONE <- function(X, D_X, phi, alpha = 1.001, beta = 1.001,
                    n.iter = 10000, thres = 1e-5){
  G <- length(X)
  CT <- length(phi)
  M <- unlist(lapply(phi, ncol))
  ct.names <- names(phi)
  phi_unmethylated <- lapply(phi, function(x){1-x})

  # Initialization
  mu <- rep(1/CT, CT)
  psi <- sapply(M, function(x){rep(1/x,x)})
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







#' cf_BLEND
#' @export
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
