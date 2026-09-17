
library(scales)
#' SelectMarkerCpGs
#' @export
SelectMarkerCpGs <- function(
    M, N, BinaryTree, ancestor.tree, leaf, sample_list, layer, upper = 81,
    maximum.idx, lowest.p = 0.9,min.read.depth = 15, min.rel.distance = 1.64,
    min.abs.distance = 0.3, file_name = "Estimation", min.IQR.target = 0.15,
    min.IQR.offtarget = 0.12, CTS.signature = T
){
  binary.tree <- BinaryTree
  Leaf.tree <- leaf
  tree_pi <- readRDS(paste0(file_name, "_tree_pi.RDS"))
  n_node <- length(Leaf.tree)
  pi.j0 <- tree_pi[1:n_node]
  pi.t <- tree_pi[-c(1:n_node)]
  n_node <- length(binary.tree)
  n_leaf <- length(Leaf.tree[[n_node]])
  onevrest <- list()
  lens = lengths(sample_list)
  samples <- sample_list
  samples_num <- lengths(samples)
  samples_num <- unname(samples_num)
  samples_num <- rep(seq_along(samples_num), times = samples_num)
  n_sample <- length(samples_num)

  leaf_all <- append_complements(leaf)
  n_col <- length(leaf_all)
  zero_idx <- which(sapply(leaf_all, function(x) identical(x, 0)))

  read_mat <- function(suffix){
    matrix(t(fread(paste0(file_name,"_",suffix,".txt"), sep = "\t"))
           ,ncol = n_col)
  }

  Loglik_full <- read_mat("ll")
  Loglik_truncate <- read_mat("ll.truncate")
  mu_full <- read_mat("mu")
  sigma_full <- read_mat("sigma")

  p_all <- calculate.p.gj(binary.tree, pi.j0, Loglik_full, layer)
  p.t <- calculate.p.t(pi.t, pi.j0, n_col, zero_idx, Loglik_full, Loglik_truncate, p_all[[2]])
  p <- p_all[[1]]
  p2 <- p_all[[1]]
  p2_1 <- p_all[[1]]
  p2_2 <- p_all[[1]]
  CalculateJoint <- function(j, current_p = NULL) {
    if (j == n_node) current_p <- 1
    if (length(binary.tree[[j]]) == 1) return()
    tmp <- current_p * (1 - p2[, j]) * p2[,binary.tree[[j]][1]] * p2[,binary.tree[[j]][2]]
    tmp1 <- current_p * (1 - p2[, j]) * p2[,binary.tree[[j]][1]]
    tmp2 <- current_p * (1 - p2[, j]) * p2[,binary.tree[[j]][2]]
    current_p <- (1 - p2[, j]) * current_p
    p2[, j] <<- tmp
    p2_1[, j] <<- tmp1
    p2_2[, j] <<- tmp2
    CalculateJoint(binary.tree[[j]][1], current_p)
    CalculateJoint(binary.tree[[j]][2], current_p)
  }
  CalculateJoint(n_node)

  ct <- names(sample_list)
  leaf_nodes <- which(sapply(binary.tree, length) == 1)
  p2[,leaf_nodes] <- 0
  trace_to_root <- function(parent_list, root = n_node) {
    paths <- vector("list", length(parent_list))
    for (i in seq_along(paths)) {
      current_node <- i
      path <- c(current_node)
      while (current_node != root) {
        current_node <- parent_list[[current_node]]
        path <- c(path, current_node)
      }
      paths[[i]] <- path
    }
    paths
  }
  paths <- trace_to_root(ancestor.tree)

  n <- nrow(mu_full)
  mu_T0_post <- matrix(0, nrow = n, ncol = n_leaf)
  sigma_T0_post <- matrix(0, nrow = n, ncol = n_leaf)

  for (i in 1:n_leaf) {
    index <- c(paths[[i]],(n_node+1):n_col)
    index[index %in% (paths[[i]] + n_node)] <- paths[[i]][paths[[i]] < n_node]
    mu_value <- mu_full[, index]
    sigma_value <- sigma_full[, index]
    weight <- cbind(p2[, paths[[i]]], p.t[, -1])
    mu_T0_post[, i] <- rowSums(mu_value * weight)
    sigma_T0_post[, i] <- rowSums(sigma_value * weight)
  }

  p_post = pnorm(mu_T0_post/sqrt(1+sigma_T0_post^2))
  alpha <- p_post*(1-p_post)/( pbv::pbvnorm(x = mu_T0_post/sqrt(1+sigma_T0_post^2), y = mu_T0_post/sqrt(1+sigma_T0_post^2),
                                            rho = sigma_T0_post^2/(1+sigma_T0_post^2)) - p_post^2 ) - 1
  beta_IQR = function(mu,sigma)
  {
    p = pnorm(mu/sqrt(1+sigma^2))
    alpha <- p*(1-p)/( pbv::pbvnorm(x = mu/sqrt(1+sigma^2), y = mu/sqrt(1+sigma^2), rho = sigma^2/(1+sigma^2)) - p^2 ) - 1
    alpha <- pmin(alpha, 150000)
    qbeta(0.75,p*alpha,(1-p)*alpha)- qbeta(0.25,p*alpha,(1-p)*alpha)
  }
  alpha = pmin(alpha,150000)
  p_post_sub = p_post[,samples_num]
  IQR = beta_IQR(mu_T0_post,sigma_T0_post)
  onevrest = list()
  for(i in 1:n)
  {
    p.tmp = p.t[i,]
    max.p = max(p.tmp)
    max.idx = which(p.tmp == max.p)[1]
    max.idx = max.idx - 1
    if(max.idx >= maximum.idx)
    {
      next
    }
    if(max.idx == 0 || max.p<lowest.p || abs((mu_full[i,n_node+max.idx] -  mu_full[i,max.idx])/
                                             sqrt(sigma_full[i,n_node+max.idx]^2 +  sigma_full[i,max.idx]^2))<min.rel.distance ||
       mean(N[i,])<min.read.depth ||any(IQR[i,Leaf.tree[[max.idx]]]>min.IQR.target)||any(IQR[i,-Leaf.tree[[max.idx]]]>min.IQR.offtarget)  )
    {
      next
    }else{
      if(abs((mu_full[i,n_node+max.idx] -  mu_full[i,max.idx])/
             sqrt(sigma_full[i,n_node+max.idx]^2 +  sigma_full[i,max.idx]^2))<min.rel.distance ||
         any(IQR[i,Leaf.tree[[max.idx]]]>min.IQR.target)||any(IQR[i,-Leaf.tree[[max.idx]]]>min.IQR.offtarget)||
         abs(mean(p_post[i,Leaf.tree[[max.idx]]]) - mean(p_post[i,-Leaf.tree[[max.idx]]]))<min.abs.distance )
      {
        next
      }
      ct.tmp = ct[max.idx]
      if(max.idx > n_leaf)
      {
        ct.tmp = max.idx
      }
      onevrest[[i]] = c(i,ct.tmp)
    }
  }

  onevrest = do.call(rbind,onevrest)
  onevrest = as.data.frame(onevrest)
  onevrest[,1] = as.numeric(onevrest[,1])
  idx_T0 = 1:n
  p2.T0 = p2[idx_T0,]
  mu.T0 = mu_full[idx_T0,1:n_node]
  sigma.T0 = sigma_full[idx_T0,1:n_node]
  length(idx_T0)
  T0_all = list()
  for(celltype in 1:(n_node-1))
  {

    T0_idx = list()
    for(i in idx_T0)
    {
      if(p.t[i,1]<lowest.p || mean(N[i,]) < min.read.depth)
      {
        next
      }
      g1 <- binary.tree[[celltype]][1]
      g2 <- binary.tree[[celltype]][2]
      mu <- mu_T0_post[i,]
      sigma <- sigma_T0_post[i,]
      posterior_p = pnorm(mu/sqrt(1+sigma^2))
      distance1 = abs((mu-mean(mu[Leaf.tree[[g1]]]))/sqrt(mean(sigma[Leaf.tree[[g1]]])^2 + sigma^2))
      distance2 = abs((mu-mean(mu[Leaf.tree[[g2]]]))/sqrt(mean(sigma[Leaf.tree[[g2]]])^2 + sigma^2))
      p_diff1 = abs(posterior_p - mean(posterior_p[Leaf.tree[[g1]]]))
      p_diff2 = abs(posterior_p - mean(posterior_p[Leaf.tree[[g2]]]))
      num = n_leaf - length(Leaf.tree[[g1]]) - length(Leaf.tree[[g2]]) - 5

      if(p2[i,celltype]>=lowest.p &&
         all(IQR[i,Leaf.tree[[g1]]] < min.IQR.target)&&
         all(IQR[i,Leaf.tree[[g2]]] < min.IQR.target)&&
         ((sum(distance1>min.rel.distance) >= num&&
         sum(p_diff1>min.abs.distance) >= num)||
         (sum(distance2>min.rel.distance) >= num&&
         sum(p_diff2>min.abs.distance) >= num)) &&
         sum(IQR[i,-c(Leaf.tree[[g1]],Leaf.tree[[g2]])]<min.IQR.offtarget) >= num)
        {
        if(sum(distance1>min.rel.distance) >= num&&
           sum(p_diff1>min.abs.distance) >= num)
        {
          T0_idx[[i]] = c(i,g1)
        }else{
          T0_idx[[i]] = c(i,g2)
        }
        }
        if(p2[i,celltype]<lowest.p &&
           celltype >= upper&&
           p2_1[i,celltype]>=lowest.p&&
           all(IQR[i,Leaf.tree[[g1]]] < min.IQR.target) &&
           sum(distance1>min.rel.distance) >= num &&
           sum(p_diff1>min.abs.distance) >= num &&
           sum(IQR[i,-c(Leaf.tree[[g1]])]<min.IQR.offtarget) >= num)
          {
            T0_idx[[i]] = c(i,g1)
          }
        if(p2[i,celltype]<lowest.p &&
           celltype >= upper&&
           p2_2[i,celltype]>=lowest.p&&
           all(IQR[i,Leaf.tree[[g2]]] < min.IQR.target) &&
           sum(distance2>min.rel.distance) >= num &&
           sum(p_diff2>min.abs.distance) >= num &&
           sum(IQR[i,-c(Leaf.tree[[g2]])]<min.IQR.offtarget) >= num)
          {
          T0_idx[[i]] = c(i,g2)
          }
      }
    T0_idx = do.call(rbind,T0_idx)
    if(is.null(T0_idx) || nrow(T0_idx)==0)
    {
      next
    }
    T0_idx = as.data.frame(T0_idx)
    T0_idx[,1] = as.numeric(T0_idx[,1])
    p_post = pnorm(mu_T0_post/sqrt(1+sigma_T0_post^2))
    p_post_sub = p_post[,samples_num]
    T0_all[[celltype]] = T0_idx
  }
  T0_all <- do.call(rbind,T0_all)
  onevrest$onevsrest <- 1
  T0_all$onevsrest <- 0
  signature <- rbind(onevrest,T0_all)
  colnames(signature) <- c("idx","celltype","onevsrest")
  signature <- cbind(signature, mu_T0_post[signature$idx,],sigma_T0_post[signature$idx,])
  if(CTS.signature == T)
  {
    out <- list()
    out[["CTS_signature"]]<- signature
    return(out)
  }
  signature_idx <- signature[,1]
  signature <- signature[order(signature_idx),]
  signature_idx <- signature_idx[order(signature_idx)]
  sigma_full <- sigma_full[signature_idx,]
  mu_full <- mu_full[signature_idx,]
  Loglik_full <- Loglik_full[signature_idx,]
  Loglik_truncate <- Loglik_truncate[signature_idx,]
  p2 <- p2[signature_idx,]
  p.t <- p.t[signature_idx,]
  p_Ind <- fread(paste0(file_name,"_Ind.txt"), sep = "\t")
  p_Ind <- as.data.frame(p_Ind)
  n <- length(signature_idx)
  p_subj_post <- matrix(rep(0,n*n_sample),ncol = n_sample)
  for(i in 1:n_sample)
  {
    j = samples_num[i]
    p_current <- p_Ind[, which(grepl(paste0("^Subj",i,"_"), names(p_Ind)))]
    p_current <- as.matrix(p_current)
    storage.mode(p_current) <- "numeric"
    index <- c(paths[[j]],(n_node+1):n_col)
    index[index %in% (paths[[j]] + n_node)] <- paths[[j]][paths[[j]] < n_node]
    weight = cbind(p2[, paths[[j]]], p.t[, -1])
    grp <- split(seq_along(index), index)
    weight <- do.call(
      cbind,
      lapply(grp, function(j) if (length(j) == 1) weight[, j] else rowSums(weight[, j, drop = FALSE]))
    )
    weight <- weight[, match(unique(index), names(grp)), drop = FALSE]
    p_current = cbind(p_current,0)
    p_subj_post[,i] = rowSums(p_current * weight)
  }
  signature_Ind <- cbind(signature[,1:3],p_subj_post)
  out <- list()
  out[["CTS_signature"]] <- signature
  out[["SubjSpec_signature"]] <- signature_Ind
  return(out)
}





