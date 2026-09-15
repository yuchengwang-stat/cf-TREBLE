#Input:
##M & U: vectors of methylated and unmethylated counts. The read depth for some samples can be 0
##sigma.grid & sigma.pi: the grid of sigma's and their prior probabilities
##mu.mu, mu.tau2, mu.pi: the prior means, variances, and mixture weights for the prior on mu
##ind.sigma.truncate: the indices of sigma.grid corresponding to the truncated prior. For example, if sigma.grid = c(0,0.1,0.2,1) and the truncated grid is c(0,0.1,0.2), then ind.sigma.truncate=c(1,2,3). It can be set to NULL, in which case ind.sigma.truncate=1:length(sigma.prior).
##Gauss.grid & Gauss.w: Gauss-Hermite points and weights. The number of points should be odd (likely 5 or 11). This should ALWAYS be specified
##is.large: T/F. Is the sample size "large"? If yes, it will use a faster (but slightly less accurate) version of the code. It defaults to T when n = length(M) >= 40

#Output:
##ll: the log-likelihood log{Pr(M, U)}
##ll.truncate: the log-likelihood log{Pr(M, U)} under the truncated prior for sigma
##mu & mu2: Posterior expectations of mu and mu^2
##p & p2: Posterior expectations of p and p^2
##sigma: Posterior expectation of sigma
##alpha: Posterior expectation of alpha conditional on sigma>0
##p0 & pg0: Posterior expectation of p conditional on sigma=0 & sigma>0, respectively
##pi.Inf.alpha: Posterior probability that alpha=Infinity (i.e., sigma=0)

PostExp_fast <- function(M, U, sigma.grid, sigma.pi, mu.mu, mu.tau2, mu.pi, ind.sigma.truncate=NULL, Gauss.grid=NULL, Gauss.w=NULL, is.large=sum(M+U>0)>=40, grad.tol=1e-6, max.iter.sub=30) {
  ind.use <- which(M+U>0); M <- M[ind.use]; U <- U[ind.use]
  sigma.pi <- sigma.pi/sum(sigma.pi); mu.pi <- mu.pi/sum(mu.pi)
  
  #Optimization + posteriors#
  if (!is.large) {
    tmp <- Create.Sigma.Mu(sigma.grid = sigma.grid, mu.mu = mu.mu, mu.tau2 = mu.tau2, sigma.pi = sigma.pi, mu.pi = mu.pi)
    Sigma <- tmp$Sigma; Mu.prior <- tmp$Mu.prior; Tau2.prior <- tmp$Tau2.prior; Prior <- tmp$Prior; rm(tmp)
    if (!is.null(ind.sigma.truncate)) {ind.sigma.truncate <- c(sapply(1:length(mu.mu),function(i){ind.sigma.truncate + (i-1)*length(sigma.grid)}))}
    opt <- My.1D.opt.sigma.mu.post(Sigma.all = Sigma, Mu.prior = Mu.prior, ind.0 = NULL, ind.g0 = NULL, Tau2.prior = Tau2.prior, M = M, U = U, add.repeat = T, tol.grad = grad.tol, max.iter = max.iter.sub)
    out <- Laplace.GQ.Post(Mu = opt$Mu, Sigma = Sigma, h = -opt$h, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Prior = Prior, M = M, U = U, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, add.repeat = T)
    out$n.iter.opt <- opt$n.iter
    return(out)
  }
  
  opt <- My.1D.opt.sigma.mu.post.UN(Sigma = sigma.grid, Mu.prior = mu.mu, Pi.prior = mu.pi, ind.0 = NULL, ind.g0 = NULL, Tau2.prior = mu.tau2, M = M, U = U, tol.grad = grad.tol, max.iter = max.iter.sub)
  if (any(is.nan(opt$h) | is.na(opt$h) | is.na(opt$f) | is.na(opt$g)) | opt$out==1) {return(PostExp_fast(M = M, U = U, sigma.grid = sigma.grid, sigma.pi = sigma.pi, mu.mu = mu.mu, mu.tau2 = mu.tau2, mu.pi = mu.pi, ind.sigma.truncate = ind.sigma.truncate, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, is.large = F, grad.tol = grad.tol))}
  if (any(opt$h>=0 | abs(opt$g)>grad.tol)) {return(PostExp_fast(M = M, U = U, sigma.grid = sigma.grid, sigma.pi = sigma.pi, mu.mu = mu.mu, mu.tau2 = mu.tau2, mu.pi = mu.pi, ind.sigma.truncate = ind.sigma.truncate, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, is.large = F, grad.tol = grad.tol))}
  ll <- opt$f - 1/2*log(-opt$h) + 1/2*log(2*pi)
  tmp <- exp(ll - max(ll))
  Exp.mu <- sum(tmp*opt$Mu*sigma.pi)/sum(tmp*sigma.pi)
  Exp.mu2 <- sum(tmp*opt$Mu^2*sigma.pi)/sum(tmp*sigma.pi)
  Exp.sigma <- sum(tmp*sigma.grid*sigma.pi)/sum(tmp*sigma.pi)
  tmp2 <- PandAlpha(theta = c(Exp.mu,Exp.sigma))
  ll.all <- max(ll)+log(sum(tmp*sigma.pi))
  ll.truncate <- ifelse(is.null(ind.sigma.truncate),ll.all, max(ll[ind.sigma.truncate]) + log(sum(exp( ll[ind.sigma.truncate] - max(ll[ind.sigma.truncate]) )*sigma.pi[ind.sigma.truncate]/sum(sigma.pi[ind.sigma.truncate]))) )
  return(list(ll=ll.all, ll.truncate=ll.truncate, mu=Exp.mu, mu2=Exp.mu2, sigma=Exp.sigma, p=tmp2[1], p2=tmp2[1]^2, p0=tmp2[1], pg0=tmp2[1], alpha=tmp2[2], pi.Inf.alpha=0, n.iter.opt=opt$n.iter))
}

Create.Sigma.Mu <- function(sigma.grid, mu.mu, mu.tau2, sigma.pi, mu.pi) {
  Sigma <- c(cbind(sigma.grid)%*%rbind(rep(1,length(mu.mu))))
  Mu.prior <- c(cbind(rep(1,length(sigma.grid)))%*%rbind(mu.mu))
  Tau2.prior <- c(cbind(rep(1,length(sigma.grid)))%*%rbind(mu.tau2))
  Prior <- c(cbind(sigma.pi)%*%rbind(mu.pi))
  return(list(Sigma=Sigma, Mu.prior=Mu.prior, Tau2.prior=Tau2.prior, Prior=Prior))
}


Laplace.GQ.Post <- function(Mu, Sigma, h, Mu.prior, Tau2.prior, Prior, M, U, ind.sigma.truncate=NULL, Gauss.grid=NULL, Gauss.w=NULL, ind.0=NULL, ind.g0=NULL, add.repeat = T) {
  if (add.repeat) {
    Sigma <- rep(Sigma,4)
    Mu.prior <- rep(Mu.prior,4) 
    Tau2.prior <- rep(Tau2.prior,4)
    ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma>0)
  }
  n <- length(Mu)/4; ind.mu <- 1:n; ind.p <- (n+1):(2*n); ind.p2 <- (2*n+1):(3*n); ind.alpha <- (3*n+1):(4*n)
  tmp <- rbind(rep(1,length(Gauss.grid)))
  Mu.grid <- cbind(sqrt(2)/sqrt(h))%*%rbind(Gauss.grid) + Mu
  tmp <- F.all.post(Mu = t(Mu.grid), Sigma = t(cbind(Sigma)%*%tmp), Mu.prior = t(cbind(Mu.prior)%*%tmp), Tau2.prior = t(cbind(Tau2.prior)%*%tmp), M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0) #GQ x Mu
  tmp <- sweep(x = t(tmp) - 1/2*log(h) + 1/2*log(2), MARGIN = 2, STATS = Gauss.grid^2 + log(Gauss.w), FUN = "+") #Mu x GQ
  C <- apply(tmp,1,max)
  tmp2 <- exp(tmp-C); r.tmp2 <- rowSums(tmp2)
  
  #Mu#
  Exp.Mu.Z = rowSums(tmp2[ind.mu,]*Mu.grid[ind.mu,])/r.tmp2[ind.mu]
  Exp.Mu2.Z = rowSums(tmp2[ind.mu,]*Mu.grid[ind.mu,]^2)/r.tmp2[ind.mu]
  post.probs <- exp(log(r.tmp2[ind.mu]) + C[ind.mu] - max(log(r.tmp2[ind.mu]) + C[ind.mu]))*Prior; post.probs <- post.probs/sum(post.probs)
  Exp.Mu <- sum(Exp.Mu.Z*post.probs); Exp.Mu2 <- sum(Exp.Mu2.Z*post.probs)
  Exp.sigma <- sum(Sigma[ind.mu]*post.probs)
  
  #p#
  tt <- post.probs*( r.tmp2[ind.p]*exp(C[ind.p] - C[ind.mu])/r.tmp2[ind.mu] )
  Exp.p <- sum(tt)
  Exp.p2 <- sum(post.probs*( r.tmp2[ind.p2]*exp(C[ind.p2] - C[ind.mu])/r.tmp2[ind.mu] ))
  ind.use <- which(!is.na(C[ind.alpha])); ind.na <- which(is.na(C[ind.alpha]))
  Exp.p0 <- ifelse(length(ind.na), sum(tt[ind.na])/sum(post.probs[ind.na]), NA)
  Exp.pg0 <- ifelse(length(ind.use), sum(tt[ind.use])/sum(post.probs[ind.use]), NA)
  
  #alpha#
  Exp.alpha <- sum(post.probs[ind.use]*( r.tmp2[ind.alpha]*exp(C[ind.alpha] - C[ind.mu])/r.tmp2[ind.mu] )[ind.use])/sum(post.probs[ind.use])
  pi.Inf.alpha <- ifelse(length(ind.na),sum(post.probs[ind.na]),NA)
  
  ##Log-likelihood##
  ll <- log(r.tmp2[ind.mu]) + C[ind.mu] + log(Prior)
  ll.all <- max(ll) + log(sum(exp(ll-max(ll))))
  if (!is.null(ind.sigma.truncate)) {
    ll.truncate <- log(r.tmp2[ind.mu[ind.sigma.truncate]]) + C[ind.mu[ind.sigma.truncate]] + log(Prior[ind.sigma.truncate]) - log(sum(Prior[ind.sigma.truncate]))
    ll.truncate <- max(ll.truncate) + log(sum(exp(ll.truncate - max(ll.truncate))))
  } else {
    ll.truncate <- ll.all
  }
  
  return(list(ll=ll.all, ll.truncate=ll.truncate, mu=Exp.Mu, mu2=Exp.Mu2, sigma=Exp.sigma, p=Exp.p, p2=Exp.p2, p0=Exp.p0, pg0=Exp.pg0, alpha=Exp.alpha, pi.Inf.alpha=pi.Inf.alpha))
}

Laplace.GQ.p <- function(Mu, Sigma, f, h, Mu.prior, Tau2.prior, ll, post.probs, M, U, Gauss.grid=NULL, Gauss.w=NULL, type=c("p","p2")) {
  type <- match.arg(type,c("p","p2"))
  tmp <- rbind(rep(1,length(Gauss.grid)))
  tmp <- F.all(Mu = cbind(sqrt(2)/sqrt(h))%*%rbind(Gauss.grid) + Mu, Sigma = cbind(Sigma)%*%tmp, Mu.prior = cbind(Mu.prior)%*%tmp, Tau2.prior = cbind(Tau2.prior)%*%tmp, M = M, U = U, ind.0 = NULL, ind.g0 = NULL, add.term = type) #Mu x GQ
  tmp <- sweep(x = tmp - 1/2*log(h) + 1/2*log(2), MARGIN = 2, STATS = Gauss.grid^2 + log(Gauss.w), FUN = "+")
  C <- apply(tmp,1,max); C.prime <- max(ll)
  return(rowSums(exp(tmp-C))/exp(ll-C.prime)*exp(C-C.prime)*post.probs)
}

Laplace.GQ.alpha <- function(Mu, Sigma, f, h, Mu.prior, Tau2.prior, ll, post.probs, M, U, Gauss.grid=NULL, Gauss.w=NULL) {
  tmp <- rbind(rep(1,length(Gauss.grid)))
  tmp <- F.all(Mu = cbind(sqrt(2)/sqrt(h))%*%rbind(Gauss.grid) + Mu, Sigma = cbind(Sigma)%*%tmp, Mu.prior = cbind(Mu.prior)%*%tmp, Tau2.prior = cbind(Tau2.prior)%*%tmp, M = M, U = U, ind.0 = NULL, ind.g0 = NULL, add.term = "alpha") #Mu x GQ
  tmp <- sweep(x = tmp - 1/2*log(h) + 1/2*log(2), MARGIN = 2, STATS = Gauss.grid^2 + log(Gauss.w), FUN = "+")
  ind.g0 <- which(Sigma>0); ind.0 <- which(Sigma==0)
  C <- apply(tmp[ind.g0,],1,max); C.prime <- max(ll[ind.g0])
  return(list(alpha=rowSums(exp(tmp[ind.g0,]-C[ind.g0]))/exp(ll[ind.g0]-C.prime)*exp(C[ind.g0]-C.prime)*post.probs[ind.g0]/sum(post.probs[ind.g0]), prob.0=ifelse(length(ind.0),sum(post.probs[ind.0]),NULL)))
}