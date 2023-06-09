library(igraph)
library(stats)
library(Matrix)
library(mc2d)

# simulation_network <- function(m, n, p=0){
#   x <- c()
#   if (p==0){
#     p <- runif(1,0,0.6) #Ideally, if p is lower than 0.5, the correspondent graph will be sparse
#   }

#   for(j in 1:m){
#     d <- rbern(n*n, p)
#     layer <- matrix(d,n,n)
#     layer <- forceSymmetric(layer, "U")
#     layer <- layer - diag(diag(layer))
#     
#     x <- append(x, as.vector(layer))
#   }
#   Y_o <- array(x, c(n, n,m)) #To be improved!!!
#   Y <- array(0,c(m,n,n))
#   for(i in 1:m){
#     Y[i,,]<-Y_o[,,i]
#   }
#   return(Y)
#   
# }

#Loss function
loss_function <- function(Y,B,h){
  #Inputs:
  # Y: Adjacency tensor w/ 0 the diagonal entries
  # B: Stochastic block model (tensor)
  # h: membership vector
  #Output:
  # Loss function
  n <- dim(Y)[3]
  B1 <- B[,h,h]
  for(i in 1:n){
    B1[,i,i] <- 0
  }
  return(sum((Y - B1)^2)) 
  
}


Algo_1_new <- function(Y, K, iter_max = 150){
  #Inputs: Y: Adjacency tensor
  #      : K: Number of communities
  
  #Outputs: g: Membership vector
  #       : B: community wise connectivity (tensor)
  
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  #Initialization
  ind <- matrix(nrow= n , ncol=m*n )
  for (j in 1:n){
    ind[j,] <- c(Y[,j,])
  }
  
  comm <- kmeans(ind, K, iter.max = 100)
  g_old <- comm$cluster
  B_old <- array(0,c(m,K,K))
  si <- loss_function(Y,B_old, g_old)
  
  B_old <- get_B(Y, g_old, K, m)
  
  #Loop until convergence criterion
  #Step 2 (loop for to determine the argmin and derive g_new):
  o = 0
  repeat{
    o <- o + 1
    g_new <- update_g(Y, B_old, g_old, K, n)
    
    #Step3: (derive B_new)
    B_new <- get_B(Y, g_new, K, m)
    
    #Step 4:
    sn <- loss_function(Y, B_new, g_new)
    #print(o)
    
    if((si > sn)&&(o<iter_max)) {
      #print(o)
      B_old <- B_new
      g_old <- g_new
    }else{
      return(g_old)
      break
    }
  }
}

get_B<- function(Y, g_vec, K, m){
  B_mat <- array(0,c(m,K,K))
  #Number of nodes in each node cluster 1 to K.
  w_vec <- sapply(1:K, function(x,g) sum(x==g), g_vec)
  #Numbr of all possible edges 
  s_mat <- tcrossprod(w_vec,w_vec)
  diag(s_mat) <- w_vec*(w_vec-1)
  
  for(k2 in 1:K){
    ind_k2 <- which(g_vec == k2)
    B_mat[, k2, k2] <- rowSums(Y[,ind_k2, ind_k2])/s_mat[k2, k2]
    for(k1 in 1:(k2-1)){
      ind_k1 <- which(g_vec == k1)
      B_mat[, k1, k2] <- rowSums(Y[,ind_k1, ind_k2])/s_mat[k1, k2]
      B_mat[, k2, k1] <- B_mat[, k1, k2] 
    }
  }
  
  return(B_mat)
}
update_g<- function(Y, B_old, g_old, K, n){
  g_new <- rep(NA, n) #vector of {1,..., K}
  for(i in 1:n){
    gm <- 10^99
    for(k in 1:K){
      gk <- sum((Y[,i,-i]-B_old[,k,g_old[-i]])^2)
      if (gk < gm){
        lk <-  k
        gm <- gk
      }
    }
    g_new[i] <- lk
  }
  return(g_new)
}
