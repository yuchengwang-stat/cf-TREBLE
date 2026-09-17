# Tree-based EM algorithm.
calculate.p.gj <- function(node.binary, pi.j, log.lik, layer) {
  log.lik.split <- log.lik
  inductive <- function(index,log.lik.split) {
    current.node <- node.binary[[index]]
    if (length(current.node) != 1) {
      left <- current.node[1]
      right <- current.node[2]
      summand <- cbind(log.lik.split[,left] + log.lik.split[,right] + log(1 - pi.j[left]) + log(1 - pi.j[right]),
                       log.lik.split[,left] + log.lik[,right] + log(1 - pi.j[left]) + log(pi.j[right]),
                       log.lik[,left] + log.lik.split[,right] + log(pi.j[left]) + log(1 - pi.j[right]),
                       log.lik[,left] + log.lik[,right] + log(pi.j[left]) + log(pi.j[right]))
      C <- apply(summand, MARGIN = 1, max)
      summand <- log(apply(exp(summand - C), MARGIN = 1, sum)) + C
      #log.lik.split[,index] <- summand
    }
    summand
  }
  
  for(k in 1:length(layer))
  {
    log.lik.split[,layer[[k]]] = mapply(inductive, layer[[k]], MoreArgs =list(log.lik.split))
  }
  

  n <- length(node.binary)
  out <- sweep((log.lik.split - log.lik)[,1:n], 2, log(1-pi.j) - log(pi.j), "+")
  out <- 1/(1+exp(out))
  list(out,log.lik.split[,n])
}
calculate.w.gj <- function(node.binary, p){
  out <- 1 - p
  recursive <- function(current.node, current.index, current.result)
  {
    if(length(current.node) == 1)
    {
      out[,current.index] <<- current.result
    }else{
      current.result.new <- current.result * out[,current.index]
      out[,current.index] <<- current.result
      left.index <- current.node[1]
      left.node <- node.binary[[left.index]]
      recursive(left.node, left.index, current.result.new)
      right.index <- current.node[2]
      right.node <- node.binary[[right.index]]
      recursive(right.node, right.index, current.result.new)
    }
    
  }
  current.index <- length(node.binary)
  current.node <- node.binary[[current.index]]
  current.result <- rep(1,nrow(out))
  recursive(current.node, current.index, current.result)
  return(out)
}
calculate.p.t = function(pi.t, pi.z, n_col, zero_idx, log.lik.full,log.lik.truncated, log.lik.root.1){
  n_1vsrest <- length(pi.t) - 1
  n_node <- length(pi.z)
  log.lik.Tk = log.lik.truncated[,(n_node+1):n_col] + 
    log.lik.truncated[,1:(n_node-1)]
  log.lik.Tk[,zero_idx - n_node] <- -Inf
  log.term.1 = log.lik.full[,n_node] + log(pi.z[n_node])
  log.term.2 = log.lik.root.1 + log(1-pi.z[n_node])
  C = pmax(log.term.1,log.term.2)
  log.lik.T0 = C + log(exp(log.term.1 - C)  + exp(log.term.2 - C))
  unnorm.log.post.T0 = log.lik.T0 + log(pi.t[1])
  unnorm.log.post.Tk = sweep(log.lik.Tk, 2, log(pi.t[-c(1)]), "+")
  unnorm.log.post.T = cbind(unnorm.log.post.T0, unnorm.log.post.Tk)
  CC = apply(unnorm.log.post.T, 1, max)
  norm.const <- log(rowSums(exp(unnorm.log.post.T - CC))) + CC
  #norm.const = CC + log(rowSums(exp(sweep(unnorm.log.post.T, 1, CC, "-"))))
  norm.log.post.T = unnorm.log.post.T - norm.const
  norm.post.T = exp(norm.log.post.T)
  return(norm.post.T)
}
calculate.pi.j <- function(p, w, p.t0){
  out <- colSums(p.t0 * p * w) / colSums(p.t0 * w)
  out
}


