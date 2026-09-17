F.all.post.UN <- function(Mu, Sigma, Mu.prior, Tau2.prior, Pi.prior=NULL, M, U, ind.0=NULL, ind.g0=NULL) {
  if (is.matrix(Mu)) {
    return(matrix(data = F.all.post.UN(Mu = c(Mu), Sigma = c(Sigma), Mu.prior = c(Mu.prior), Tau2.prior = c(Tau2.prior), Pi.prior = NULL, M = M, U = U, ind.0 = NULL, ind.g0 = NULL), ncol = ncol(Mu)))
  }
  if (is.null(ind.0)) {ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma>0)}
  if (!is.null(Pi.prior)) {
    out <- Mixture.Gauss(Mu = Mu, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Pi.prior = Pi.prior, deriv = 0)
  } else {
    out <- dnorm(x = Mu, mean = Mu.prior, sd = sqrt(Tau2.prior), log = T)
  }
  if (length(ind.0)) {
    p <- pnorm(Mu[ind.0])
    out[ind.0] <- out[ind.0] + sum(M)*log(p) + sum(U)*log(1-p)
  }
  if (length(ind.g0)) {
    tmp <- PandAlpha(theta = cbind(Mu[ind.g0],Sigma[ind.g0]))
    p <- tmp[,1]; alpha <- tmp[,2]
    a <- p*alpha; b <- (1-p)*alpha
    n <- length(M)
    out[ind.g0] <- out[ind.g0] + rowSums(lbeta(cbind(a)%*%rbind(rep(1,n)) + cbind(rep(1,length(a)))%*%rbind(M), cbind(b)%*%rbind(rep(1,n)) + cbind(rep(1,length(a)))%*%rbind(U))) - n*lbeta(a, b)
  }
  return(out)
}

Grad.F.all.post.UN <- function(Mu, Sigma, Mu.prior, Tau2.prior, Pi.prior=NULL, M, U, ind.0, ind.g0) {
  if (!is.null(Pi.prior)) {
    out <- Mixture.Gauss(Mu = Mu, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Pi.prior = Pi.prior, deriv = 1)
  } else {
    out <- -(Mu - Mu.prior)/Tau2.prior
  }
  
  if (length(ind.g0)) {
    tmp <- Get.Jacob(theta = cbind(Mu[ind.g0],Sigma[ind.g0]), return.PandAlpha = T)
    J <- rbind(tmp$J); p <- rbind(tmp$p)[,1]; alpha <- rbind(tmp$p)[,2]
    a <- p*alpha; b <- (1-p)*alpha; N <- length(M)
    tmp1 <- rowSums(digamma(cbind(alpha)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M+U)))
    tmp2 <- N*digamma(alpha)
    dhda <- rowSums(digamma(cbind(a)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M))) - tmp1 - N*digamma(a) + tmp2
    dhdb <- rowSums(digamma(cbind(b)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(U))) - tmp1 - N*digamma(b) + tmp2
    grad.h <- cbind(dhda*alpha - dhdb*alpha, dhda*p + dhdb*(1-p))
    
    out[ind.g0] <- out[ind.g0] + rowSums(J*grad.h)
  }
  if (length(ind.0)) {
    p <- pnorm(Mu[ind.0])
    out[ind.0] <- out[ind.0] + dnorm(Mu[ind.0])*(sum(M)/p - sum(U)/(1-p))
  }
  return(out)
}

Hess.F.all.post.UN <- function(Mu, Sigma, Mu.prior, Tau2.prior, Pi.prior=NULL, M, U, ind.0, ind.g0) {
  if (!is.null(Pi.prior)) {
    out <- Mixture.Gauss(Mu = Mu, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Pi.prior = Pi.prior, deriv = 2)
  } else {
    out <- -1/Tau2.prior
  }
  
  if (length(ind.g0)) {
    tmp <- Grad.J(theta = cbind(Mu[ind.g0],Sigma[ind.g0]))
    J <- tmp$J; p <- tmp$p; alpha <- tmp$alpha; dJ.dmu <- tmp$dJ.dmu
    a <- p*alpha; b <- (1-p)*alpha
    N <- length(M)
    grad.h <- Grad.h(theta = cbind(p,alpha), M = M, U = U)
    
    tmp.1 <- rowSums(trigamma(cbind(a)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M)))
    tmp.2 <- rowSums(trigamma(cbind(b)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(U)))
    tmp.3 <- trigamma(a); tmp.4 <- trigamma(b)
    H.11 <- alpha^2*(tmp.1 + tmp.2) - N*alpha^2*tmp.3 - N*alpha^2*tmp.4
    H.22 <- p^2*tmp.1 + (1-p)^2*tmp.2 - rowSums(trigamma(cbind(alpha)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M+U))) - N*p^2*tmp.3 - N*(1-p)^2*tmp.4 + N*trigamma(alpha)
    H.12 <- rowSums(digamma(cbind(a)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M))) - rowSums(digamma(cbind(b)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(U))) - N*digamma(a) + N*digamma(b)
    H.12 <- H.12 + alpha*p*tmp.1 - alpha*(1-p)*tmp.2 - N*alpha*p*tmp.3 + N*alpha*(1-p)*tmp.4
    
    out[ind.g0] <- out[ind.g0] + J[,1]^2*H.11 + J[,2]^2*H.22 + 2*J[,1]*J[,2]*H.12 + rowSums(dJ.dmu*grad.h)
  }
  if (length(ind.0)) {
    p <- pnorm(Mu[ind.0])
    M <- sum(M); U <- sum(U)
    out[ind.0] <- out[ind.0] + phi.prime(x = Mu[ind.0])*(M/p - U/(1-p)) - dnorm(Mu[ind.0])^2*(M/p^2 + U/(1-p)^2)
  }
  return(out)
}

My.1D.opt.sigma.mu.post.UN <- function(Sigma, Mu.prior, ind.0=NULL, ind.g0=NULL, Tau2.prior, Pi.prior, M, U, max.iter=1e4, tol.grad=1e-6, alpha=0.5, beta=0.5) {
  if (is.null(ind.0) | is.null(ind.g0)) {ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma>0)}
  Mu <- rep(qnorm((sum(M)+alpha)/(sum(M+U)+alpha+beta)), length(Sigma))
  for (i in 1:max.iter) {
    g <- Grad.F.all.post.UN(Mu = Mu, Sigma = Sigma, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Pi.prior = Pi.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0)
    h <- Hess.F.all.post.UN(Mu = Mu, Sigma = Sigma, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, Pi.prior = Pi.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0)
    if (any(is.na(g) | is.na(h))) {return(list(Mu=Mu, g=g, h=h, f=F.all.post.UN(Mu = Mu, Sigma = Sigma, Mu.prior = Mu.prior, Pi.prior = Pi.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0), out=0, n.iter=i))}
    if (max(abs(g),na.rm=T)<tol.grad) {
      return(list(Mu=Mu, g=g, h=h, f=F.all.post.UN(Mu = Mu, Sigma = Sigma, Mu.prior = Mu.prior, Pi.prior = Pi.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0), out=0, n.iter=i))
    }
    Mu <- Mu - g/h
  }
  return(list(Mu=Mu, g=g, h=h, f=F.all.post.UN(Mu = Mu, Sigma = Sigma, Mu.prior = Mu.prior, Pi.prior = Pi.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0), out=1, n.iter=i))
}

Mixture.Gauss <- function(Mu, Mu.prior, Tau2.prior, Pi.prior, deriv=c(0,1,2)) {
  tmp <- dnorm(x = cbind(Mu)%*%rbind(rep(1,length(Mu.prior))), mean = cbind(rep(1,length(Mu)))%*%rbind(Mu.prior), sd = cbind(rep(1,length(Mu)))%*%rbind(sqrt(Tau2.prior)), log = T)
  C <- apply(tmp,1,max)
  tmp <- exp(tmp - C)
  denom <- c(tmp%*%Pi.prior)
  if (deriv==0) {return(C + log(denom))}
  f.prime <- -( cbind(Mu)%*%rbind(rep(1,length(Mu.prior))) - cbind(rep(1,length(Mu)))%*%rbind(Mu.prior) )/(cbind(rep(1,length(Mu)))%*%rbind(Tau2.prior))
  deriv.1 <- c((tmp*f.prime)%*%Pi.prior)/denom
  if (deriv==1) {return(deriv.1)}
  f.prime2 <- -(cbind(rep(1,length(Mu)))%*%rbind(1/Tau2.prior))
  return(c((tmp*(f.prime^2 + f.prime2))%*%Pi.prior)/denom - deriv.1^2)
}