PostExp_IndSig <- function(M, U, sigma.grid, sigma.pi, mu.mu, mu.tau2, mu.pi, Gauss.grid=NULL, Gauss.w=NULL, is.large=sum(M+U>0)>=20, grad.tol=1e-6, max.iter.sub=30, return.all=F, ind.use=NULL, length.M=NULL) {
  out.post <- rep(NA,ifelse(is.null(length.M),length(M),length.M))
  if (is.null(ind.use)) {
    ind.use <- which(M+U>0)
    M <- M[ind.use]; U <- U[ind.use]
  }
  if (!length(ind.use)) {return(out.post)}
  sigma.pi <- sigma.pi/sum(sigma.pi); mu.pi <- mu.pi/sum(mu.pi)
  
  #Optimization + posteriors#
  if (!is.large) {
    tmp <- Create.Sigma.Mu(sigma.grid = sigma.grid, mu.mu = mu.mu, mu.tau2 = mu.tau2, sigma.pi = sigma.pi, mu.pi = mu.pi)
    Sigma <- tmp$Sigma; Mu.prior <- tmp$Mu.prior; Tau2.prior <- tmp$Tau2.prior; Prior <- tmp$Prior; rm(tmp)
    ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma > 0)
    opt <- My.1D.opt.sigma.mu.post.UN(Sigma = Sigma, Mu.prior = Mu.prior, ind.0 = ind.0, ind.g0 = ind.g0, Tau2.prior = Tau2.prior, M = M, U = U, Pi.prior = NULL, max.iter = max.iter.sub, tol.grad = grad.tol)
    out.post[ind.use] <- Laplace.GQ.Post_IndSig(Mu = opt$Mu, Sigma = Sigma, h = -opt$h, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Prior = Prior, M = M, U = U, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, ind.0 = ind.0, ind.g0 = ind.g0)
    return(out.post)
  }
  
  opt <- My.1D.opt.sigma.mu.post.UN(Sigma = sigma.grid, Mu.prior = mu.mu, Pi.prior = mu.pi, ind.0 = NULL, ind.g0 = NULL, Tau2.prior = mu.tau2, M = M, U = U, tol.grad = grad.tol, max.iter = max.iter.sub)
  if (any(is.nan(opt$h) | is.na(opt$h) | is.na(opt$f) | is.na(opt$g)) | opt$out==1) {return(PostExp_IndSig(M = M, U = U, sigma.grid = sigma.grid, sigma.pi = sigma.pi, mu.mu = mu.mu, mu.tau2 = mu.tau2, mu.pi = mu.pi, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, is.large = F, grad.tol = grad.tol, ind.use = ind.use, length.M = length(out.post)))}
  if (any(opt$h>=0 | abs(opt$g)>grad.tol)) {return(PostExp_IndSig(M = M, U = U, sigma.grid = sigma.grid, sigma.pi = sigma.pi, mu.mu = mu.mu, mu.tau2 = mu.tau2, mu.pi = mu.pi, Gauss.grid = Gauss.grid, Gauss.w = Gauss.w, is.large = F, grad.tol = grad.tol, ind.use = ind.use, length.M = length(out.post)))}
  if (return.all) {return(Laplace.Post_IndSig(Mu = opt$Mu, Sigma = sigma.grid, f = opt$f, h = -opt$h, Prior = sigma.pi, M = M, U = U, return.all=return.all))}
  out.post[ind.use] <- Laplace.Post_IndSig(Mu = opt$Mu, Sigma = sigma.grid, f = opt$f, h = -opt$h, Prior = sigma.pi, M = M, U = U, return.all=F)
  return(out.post)
}

Laplace.GQ.Post_IndSig <- function(Mu, Sigma, h, Mu.prior, Tau2.prior, Prior, M, U, Gauss.grid, Gauss.w, ind.0, ind.g0) {
  tmp.vec <- rbind(rep(1,length(Gauss.grid)))
  Mu.grid <- cbind(sqrt(2)/sqrt(h))%*%rbind(Gauss.grid) + Mu
  tmp <- F.all.post.UN(Mu = t(Mu.grid), Sigma = t(cbind(Sigma)%*%tmp.vec), Mu.prior = t(cbind(Mu.prior)%*%tmp.vec), Tau2.prior = t(cbind(Tau2.prior)%*%tmp.vec), M = M, U = U) #GQ x Mu
  tmp <- sweep(x = t(tmp) - 1/2*log(h) + 1/2*log(2), MARGIN = 2, STATS = Gauss.grid^2 + log(Gauss.w), FUN = "+") #Mu x GQ
  tmp.num <- c(tmp) + postexp.p(Mu = c(Mu.grid), Sigma = c(cbind(Sigma)%*%tmp.vec), M = M, U = U, use.log = T)
  C <- apply(tmp,1,max)
  
  ll.denom <- log(rowSums(exp(tmp-C))) + C  #log denominator
  post.probs <- exp(ll.denom - max(ll.denom))*Prior; post.probs <- post.probs/sum(post.probs)
  return(apply(cbind(tmp.num),2,function(x){
    x <- matrix(data = x, nrow = nrow(tmp))
    C <- apply(x,1,max)
    ll <- log(rowSums(exp(x-C))) + C
    return(sum(exp(ll - ll.denom + log(post.probs))))
  }))  
}

Laplace.Post_IndSig <- function(Mu, Sigma, f, h, Prior, M, U, return.all=F) {
  ll.denom <- f - 1/2*log(h) + 1/2*log(2*pi)
  post.probs <- exp(ll.denom - max(ll.denom))*Prior; post.probs <- post.probs/sum(post.probs)
  if (!return.all) {return(c(t(postexp.p(Mu = Mu, Sigma = Sigma, M = M, U = U, use.log = F))%*%post.probs))}
  return(list(post.probs=post.probs,Post.Exp.all=t(postexp.p(Mu = Mu, Sigma = Sigma, M = M, U = U, use.log = F)),Post.Exp=c(t(postexp.p(Mu = Mu, Sigma = Sigma, M = M, U = U, use.log = F))%*%post.probs)))
}

postexp.p <- function(Mu, Sigma, M, U, use.log=T) {
  Mu <- c(Mu); Sigma <- c(Sigma)
  ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma>0); n <- length(M)
  out <- matrix(NA,nrow=length(Mu),ncol=n)
  if (length(ind.0)==1 | (length(ind.0)>1 & n==1)) {out[ind.0,] <- pnorm(Mu[ind.0])}
  if (length(ind.0)>1 & n>1) {out[ind.0,] <- cbind(pnorm(Mu[ind.0]))%*%rbind(rep(1,n))}
  if (length(ind.g0)) {
    tmp <- PandAlpha(theta = cbind(Mu[ind.g0],Sigma[ind.g0]))
    if (n==1) {
      out[ind.g0,] <- (M + tmp[,2]*tmp[,1])/(M + U + tmp[,2])
    } else {
      out[ind.g0,] <- ( cbind(rep(1,length(ind.g0)))%*%rbind(M) + cbind(tmp[,2]*tmp[,1])%*%rbind(rep(1,n)) )/( cbind(rep(1,length(ind.g0)))%*%rbind(M+U) + cbind(tmp[,2])%*%rbind(rep(1,n)) )
    }
  }
  if (n==1) {out <- c(out)}
  if (use.log) {out <- log(out)}
  return(out)
}
