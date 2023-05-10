library(rTensor)
library(ggplot2)
library(plotly)
library(Matrix)
library(rMultiNet)


norm_vec <- function(x) sqrt(sum(x^2))
reg_vec <- function(x,delta) min(delta,norm_vec(x))/norm_vec(x)*x

PowerIteration<- function(tnsr, ranks=NULL, type="TWIST", U_0_list, delta1=1000, delta2=1000, max_iter = 25, tol = 1e-05){
  stopifnot(is(tnsr, "Tensor"))
  if (is.null(ranks))
    stop("ranks must be specified")
  if (sum(ranks > tnsr@modes) != 0)
    stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  if(type == "TWIST")
  {
    num_modes <- tnsr@num_modes
    U_list <- U_0_list
    tnsr_norm <- fnorm(tnsr)
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- rep(0, max_iter)
    CHECK_CONV <- function(Z, U_list) {
      est <- ttl(Z, U_list, ms = 1:num_modes)
      curr_resid <- fnorm(tnsr - est)
      fnorm_resid[curr_iter] <<- curr_resid
      if (curr_iter == 1)
        return(FALSE)
      if (abs(curr_resid - fnorm_resid[curr_iter - 1])/tnsr_norm <
          tol)
        return(TRUE)
      else {
        return(FALSE)
      }
    }
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
    while ((curr_iter < max_iter) && (!converged)) {
      cat("iteration", curr_iter, "\n")
      setTxtProgressBar(pb, curr_iter)
      modes <- tnsr@modes
      modes_seq <- 1:num_modes

      ##Regularization
      U_list_reg = U_list
      for(m in modes_seq)
      {
        if(m == 1 | m == 2)
        {
          U_list_reg[[m]] = t(apply(U_list_reg[[m]],1,reg_vec, delta=delta1))
        }
        if(m == 3)
        {
          U_list_reg[[m]] = t(apply(U_list_reg[[m]],1,reg_vec, delta=delta2))
        }
      }

      ##Iterate
      for (m in modes_seq) {
        X <- ttl(tnsr, lapply(U_list_reg[-m], t), ms = modes_seq[-m])
        U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
      }

      Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)

      if (CHECK_CONV(Z, U_list)) {
        converged <- TRUE
        setTxtProgressBar(pb, max_iter)
      }
      else {
        curr_iter <- curr_iter + 1
      }
    }
    close(pb)
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1)/tnsr_norm)) *
      100
    est <- ttl(Z, U_list, ms = 1:num_modes)
    invisible(list(Z = Z, U = U_list, conv = converged, est = est,
                   norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid,
                                                                   1), all_resids = fnorm_resid))
    network_embedding <- U_list[[3]]
    node_embedding <- U_list[[1]]
    return(list(Z, network_embedding, node_embedding))
  }
  if(type == "TUCKER")
  {
    decomp=tucker(arrT,ranks,max_iter = 10000,tol=1e-05)
    nodes_embedding_Our=decomp[["U"]][[1]] #nodes' embedding
    network_embedding=decomp[["U"]][[3]] #layers' embedding
    Z = decomp[["Z"]]
    return(list(Z, network_embedding, node_embedding))
  }

}

layer_comm_2 <- function(W, m){
  l <- dim(W)[1]
  g <- kmeans(W, m)$cluster
  L <- 1:m
  Z <- array(0, c(l, m))
  for(i in 1:l){
    Z[i,] <- c(L== g[i])
  }
  
  return(list(g, Z))
}

nodes_comm_2 <- function(A, g_lay, K){
  A <- A@data
  n <- dim(A)[2]
  m <- length(K)
  Z <- list()
  g_nodes <- list()
  for(i in 1:m){
    l <- which(g_lay == i)
    Al <- matrix(A[l,,],nrow = n)
    g_nodes[[i]] <- kmeans(Al, K[i])$cluster
    Zi <- array(0, c(n, m))
    L <- 1:m
    for(j in 1:n){
      Zi[j,] <- c(L== g_nodes[[i]][j])
    }
    Z[[i]] <- Zi
  }
  return(list(g_nodes,Z))
  
}
# #Number of inidividuals
# N <- 64
# #Communities
# K <- c(2,2)
# #Class of layer
# L <- c(4,4)
# M <- 2
# 
# netwe <- simul_final(N,K,L,M)
# A <- array(0,c(N,N,8))
# for(i in 1:8){
#   A[,,i] <- netwe[[1]][i,,]
# }
# tnsr <- as.tensor(A)
# r <- c(8,8,2)
# U_init <- InitializationMMSBM(tnsr, ranks = r )
# f <- PowerIteration(tnsr, ranks = r, type="TWIST", U_init, delta1=1000, delta2=1000, max_iter = 25, tol = 1e-05)
# print(dim(f[[3]]))
# # layes <- layer_comm_2(f[[2]], M)
# # nod <- nodes_comm_2(A, layes[[1]], K)



