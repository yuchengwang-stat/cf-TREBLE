F.all.post <- function(Mu, Mu.p, Mu.p2, Mu.alpha, Sigma, Mu.prior, Tau2.prior, M, U, ind.0=NULL, ind.g0=NULL, n.col=NA) {
  if (is.vector(Mu)) {
    if (is.null(ind.0)) {ind.0 <- which(Sigma==0); ind.g0 <- which(Sigma>0)}
    out <- -1/2*log(Tau2.prior*2*pi) - 1/2/Tau2.prior*(Mu - Mu.prior)^2
    n.0 <- length(ind.0)/4; n.g0 <- length(ind.g0)/4
    if (n.0) {
      ind.0.p <- ind.0[(n.0+1):(2*n.0)]; ind.0.p2 <- ind.0[(2*n.0+1):(3*n.0)]; ind.0.alpha <- ind.0[(3*n.0+1):(4*n.0)]
      p <- pnorm(Mu[ind.0])
      out[ind.0] <- out[ind.0] + sum(M)*log(p) + sum(U)*log(1-p)
      out[ind.0.p] <- out[ind.0.p] + log(p[(n.0+1):(2*n.0)])
      out[ind.0.p2] <- out[ind.0.p2] + 2*log(p[(2*n.0+1):(3*n.0)])
      out[ind.0.alpha] <- NA
    }
    if (n.g0) {
      ind.g0.p <- ind.g0[(n.g0+1):(2*n.g0)]; ind.g0.p2 <- ind.g0[(2*n.g0+1):(3*n.g0)]; ind.g0.alpha <- ind.g0[(3*n.g0+1):(4*n.g0)]
      tmp <- PandAlpha(theta = cbind(Mu[ind.g0],Sigma[ind.g0]))
      p <- tmp[,1]; alpha <- tmp[,2]
      a <- p*alpha; b <- (1-p)*alpha
      n <- length(M)
      out[ind.g0] <- out[ind.g0] + rowSums(lbeta(cbind(a)%*%rbind(rep(1,n)) + cbind(rep(1,length(a)))%*%rbind(M), cbind(b)%*%rbind(rep(1,n)) + cbind(rep(1,length(a)))%*%rbind(U))) - n*lbeta(a, b)
      out[ind.g0.p] <- out[ind.g0.p] + log(p[(n.g0+1):(2*n.g0)])
      out[ind.g0.p2] <- out[ind.g0.p2] + 2*log(p[(2*n.g0+1):(3*n.g0)])
      out[ind.g0.alpha] <- out[ind.g0.alpha] + log(alpha[(3*n.g0+1):(4*n.g0)])
    }
    if (is.na(n.col)) {return(out)}
    return(matrix(data = out, ncol = n.col))
  }
  return(F.all.post(Mu = c(Mu), Sigma = c(Sigma), Mu.prior = c(Mu.prior), Tau2.prior = c(Tau2.prior), M = M, U = U, ind.0 = which(c(Sigma)==0), ind.g0 = which(c(Sigma)>0), n.col = ncol(Mu)))
}

Grad.F.all.post <- function(Mu, Sigma, Mu.prior, Tau2.prior, M, U, ind.0, ind.g0) {
  out <- -(Mu - Mu.prior)/Tau2.prior
  n.0 <- length(ind.0)/4; n.g0 <- length(ind.g0)/4
  
  if (n.g0) {
    ind.g0.p <- ind.g0[(n.g0+1):(2*n.g0)]; ind.g0.p2 <- ind.g0[(2*n.g0+1):(3*n.g0)]; ind.g0.alpha <- ind.g0[(3*n.g0+1):(4*n.g0)]
    tmp <- Get.Jacob(theta = cbind(Mu[ind.g0],Sigma[ind.g0]), return.PandAlpha = T)
    J <- rbind(tmp$J); p <- rbind(tmp$p)[,1]; alpha <- rbind(tmp$p)[,2]
    a <- p*alpha; b <- (1-p)*alpha; N <- length(M)
    tmp1 <- rowSums(digamma(cbind(alpha)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M+U)))
    tmp2 <- N*digamma(alpha)
    dhda <- rowSums(digamma(cbind(a)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M))) - tmp1 - N*digamma(a) + tmp2
    dhdb <- rowSums(digamma(cbind(b)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(U))) - tmp1 - N*digamma(b) + tmp2
    grad.h <- cbind(dhda*alpha - dhdb*alpha, dhda*p + dhdb*(1-p))
    
    out[ind.g0] <- out[ind.g0] + rowSums(J*grad.h)
    out[ind.g0.p] <- out[ind.g0.p] + J[(n.g0+1):(2*n.g0),1]/p[(n.g0+1):(2*n.g0)]
    out[ind.g0.p2] <- out[ind.g0.p2] + 2*J[(2*n.g0+1):(3*n.g0),1]/p[(2*n.g0+1):(3*n.g0)]
    out[ind.g0.alpha] <- out[ind.g0.alpha] + J[(3*n.g0+1):(4*n.g0),2]/alpha[(3*n.g0+1):(4*n.g0)]
  }
  if (n.0) {
    ind.0.p <- ind.0[(n.0+1):(2*n.0)]; ind.0.p2 <- ind.0[(2*n.0+1):(3*n.0)]; ind.0.alpha <- ind.0[(3*n.0+1):(4*n.0)]
    p <- pnorm(Mu[ind.0])
    out[ind.0] <- out[ind.0] + dnorm(Mu[ind.0])*(sum(M)/p - sum(U)/(1-p))
    out[ind.0.p] <- out[ind.0.p] + dnorm(Mu[ind.0.p])/p[(n.0+1):(2*n.0)]
    out[ind.0.p2] <- out[ind.0.p2] + 2*dnorm(Mu[ind.0.p2])/p[(2*n.0+1):(3*n.0)]
    out[ind.0.alpha] <- NA
  }
  return(out)
}

Hess.F.all.post <- function(Mu, Sigma, Mu.prior, Tau2.prior, M, U, ind.0, ind.g0) {
  out <- -1/Tau2.prior
  n.0 <- length(ind.0)/4; n.g0 <- length(ind.g0)/4
  
  if (n.g0) {
    ind.g0.p <- ind.g0[(n.g0+1):(2*n.g0)]; ind.g0.p2 <- ind.g0[(2*n.g0+1):(3*n.g0)]; ind.g0.alpha <- ind.g0[(3*n.g0+1):(4*n.g0)]
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
    out[ind.g0.p] <- out[ind.g0.p] + (dJ.dmu[,1]/p - J[,1]^2/p^2)[(n.g0+1):(2*n.g0)]
    out[ind.g0.p2] <- out[ind.g0.p2] + 2*(dJ.dmu[,1]/p - J[,1]^2/p^2)[(2*n.g0+1):(3*n.g0)]
    out[ind.g0.alpha] <- out[ind.g0.alpha] + (dJ.dmu[,2]/alpha - J[,2]^2/alpha^2)[(3*n.g0+1):(4*n.g0)]
  }
  if (n.0) {
    ind.0.p <- ind.0[(n.0+1):(2*n.0)]; ind.0.p2 <- ind.0[(2*n.0+1):(3*n.0)]; ind.0.alpha <- ind.0[(3*n.0+1):(4*n.0)]
    p <- pnorm(Mu[ind.0])
    M <- sum(M); U <- sum(U)
    out[ind.0] <- out[ind.0] + phi.prime(x = Mu[ind.0])*(M/p - U/(1-p)) - dnorm(Mu[ind.0])^2*(M/p^2 + U/(1-p)^2)
    foo <- phi.prime(Mu[ind.0])/p - (dnorm(Mu[ind.0])/p)^2
    out[ind.0.p] <- out[ind.0.p] + foo[(n.0+1):(2*n.0)]
    out[ind.0.p2] <- out[ind.0.p2] + 2*foo[(2*n.0+1):(3*n.0)]
    out[ind.0.alpha] <- NA
  }
  return(out)
}

My.1D.opt.sigma.mu.post <- function(Sigma.all, Mu.prior, ind.0, ind.g0, Tau2.prior, M, U, max.iter=1e4, tol.grad=1e-6, alpha=0.5, beta=0.5, add.repeat=F) {
  if (add.repeat) {
    Sigma.all <- rep(Sigma.all,4)
    Mu.prior <- rep(Mu.prior,4) 
    Tau2.prior <- rep(Tau2.prior,4)
    ind.0 <- which(Sigma.all==0); ind.g0 <- which(Sigma.all>0)
   }
  Mu <- rep(qnorm((sum(M)+alpha)/(sum(M+U)+alpha+beta)), length(Sigma.all))
  for (i in 1:max.iter) {
    g <- Grad.F.all.post(Mu = Mu, Sigma = Sigma.all, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0)
    h <- Hess.F.all.post(Mu = Mu, Sigma = Sigma.all, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0)
    if (max(abs(g),na.rm=T)<tol.grad) {
      return(list(Mu=Mu, g=g, h=h, f=F.all.post(Mu = Mu, Sigma = Sigma.all, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0), out=0, n.iter=i))
    }
    Mu <- Mu - g/h
  }
  return(list(Mu=Mu, g=g, h=h, f=F.all.post(Mu = Mu, Sigma = Sigma.all, Mu.prior = Mu.prior, Tau2.prior = Tau2.prior, M = M, U = U, ind.0 = ind.0, ind.g0 = ind.g0), out=1, n.iter=i))
}

phi.prime <- function(x) {
  return(-x*dnorm(x))
}

phi.prime2 <- function(x, mean=0, sd=1) {
  return(dnorm(x = x, mean = mean, sd = sd)*(mean-x)/sd^2)
}

Get.Jacob <- function(theta, return.PandAlpha=F) {
  if (is.vector(theta)) {
    mu <- theta[1]; sigma <- theta[2]
    p <- pnorm(mu/sqrt(1+sigma^2))
    q <- pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2))
    J <- rep(NA,2)
    J[1] <- 1/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))
    dqdmu <- 2*dnorm(x = mu, mean = 0, sd = sqrt(1+sigma^2))*pnorm(mu/sqrt( (1+sigma^2)*(1+2*sigma^2) ))
    J[2] <- (J[1] - 2*p*J[1])/(q-p^2) - p*(1-p)/(q-p^2)^2*(dqdmu - 2*p*J[1])
    if (!return.PandAlpha) {return(J)}
    return(list(J=J,p=c(p,p*(1-p)/(q-p^2) - 1)))
  }
  mu <- theta[,1]; sigma <- theta[,2]
  p <- pnorm(mu/sqrt(1+sigma^2))
  q <- pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2))
  J <- matrix(NA,nrow=length(mu),ncol=2)
  J[,1] <- 1/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))
  dqdmu <- 2*dnorm(x = mu, mean = 0, sd = sqrt(1+sigma^2))*pnorm(mu/sqrt( (1+sigma^2)*(1+2*sigma^2) ))
  J[,2] <- (J[,1] - 2*p*J[,1])/(q-p^2) - p*(1-p)/(q-p^2)^2*(dqdmu - 2*p*J[,1])  
  if (!return.PandAlpha) {return(J)}
  return(list(J=J,p=cbind(p,p*(1-p)/(q-p^2) - 1)))
}

Grad.h <- function(theta, M, U) {
  p <- theta[,1]; alpha <- theta[,2]
  a <- p*alpha; b <- (1-p)*alpha; N <- length(M)
  tmp1 <- rowSums(digamma(cbind(alpha)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M+U)))
  tmp2 <- N*digamma(alpha)
  dhda <- rowSums(digamma(cbind(a)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(M))) - tmp1 - N*digamma(a) + tmp2
  dhdb <- rowSums(digamma(cbind(b)%*%rbind(rep(1,N)) + cbind(rep(1,length(a)))%*%rbind(U))) - tmp1 - N*digamma(b) + tmp2
  return(cbind(dhda*alpha - dhdb*alpha, dhda*p + dhdb*(1-p)))
}


Grad.J <- function(theta) {
  mu <- theta[,1]; sigma <- theta[,2]; v <- sigma^2
  p <- pnorm(mu/sqrt(1+sigma^2))
  q <- pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2))
  J <- Get.Jacob(theta = theta, return.PandAlpha = F)
  dqdmu <- 2*dnorm(x = mu, mean = 0, sd = sqrt(1+sigma^2))*pnorm(mu/sqrt( (1+sigma^2)*(1+2*sigma^2) ))
  dq2dmu2 <- 2*( phi.prime2(x = mu, mean = 0, sd = sqrt(1+v))*pnorm(mu/sqrt( (1+sigma^2)*(1+2*sigma^2) )) + dnorm(x = mu, mean = 0, sd = sqrt(1+sigma^2))/sqrt( (1+sigma^2)*(1+2*sigma^2) )*dnorm(mu/sqrt( (1+sigma^2)*(1+2*sigma^2) )) )
  
  #dJ/dmu#
  dJ.dmu <- matrix(NA,ncol=2,nrow=nrow(J))
  dJ.dmu[,1] <- 1/(1+v)*phi.prime(mu/sqrt(1+v))
  dJ.dmu[,2] <- (dJ.dmu[,1] - 2*p*dJ.dmu[,1] - 2*J[,1]^2)/(q - p^2)
  dJ.dmu[,2] <- dJ.dmu[,2] - (J[,1] - 2*p*J[,1])*(dqdmu - 2*p*J[,1])/(q-p^2)^2
  dJ.dmu[,2] <- dJ.dmu[,2] - (dqdmu - 2*p*J[,1])*( (J[,1] - 2*p*J[,1])/(q-p^2)^2 - 2*(p-p^2)/(q-p^2)^3*(dqdmu - 2*p*J[,1]) )
  dJ.dmu[,2] <- dJ.dmu[,2] - p*(1-p)/(q-p^2)^2*( dq2dmu2 - 2*J[,1]^2 - 2*p*dJ.dmu[,1] )
  
  return(list(p=p, alpha=p*(1-p)/(q-p^2)-1, J=J, dJ.dmu=dJ.dmu))
}

PandAlpha <- function(theta) {
  if (is.vector(theta)) {
    mu <- theta[1]; sigma <- theta[2]
    if (sigma==0) {return(c(pnorm(mu),Inf))}
    p <- pnorm(mu/sqrt(1+sigma^2))
    alpha <- p*(1-p)/( pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2)) - p^2 ) - 1
    return(c(p,alpha))
  }
  mu <- theta[,1]; sigma <- theta[,2]
  p <- pnorm(mu/sqrt(1+sigma^2))
  alpha <- p*(1-p)/( pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2)) - p^2 ) - 1
  return(cbind(p,alpha))
}