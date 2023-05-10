library(igraph)
library(stats)
library(Matrix)
library(mc2d)
library(RSpectra)

msqrt <- function(M){
  #Inputs
  # M: a (square and pos. def.) matrix
  E <- eigen(M)
  d <- (E$values)^(-1/2)
  P <- E$vectors
  
  return(tcrossprod(P %*% diag(d),P))
  
}

projection_0m <- function(X){
  #Inputs:
  # X: a matrix
  Y <- msqrt(crossprod(X))
  return( X %*% Y)
  
}
projection_rankm <- function(A, k){ #LOW-RANK projection of matrices
  #Inputs: 
  # A: array (or matrix) of rank n > k
  # k: rank required
  n <- dim(A)[1]
  S <- svds(A,k, nu = k, nv = k)
  U <- S$u
  sigma <- S$d
  D <- diag(sigma)
  V <- S$v
  return(tcrossprod(U %*% D, V))
  
}

projection_rankt <- function(A, K){#Low-rank projection of tensors
  #Inputs:
  # A: tensor, where A[i,,] is a matrix of rank ni > ki
  # k: vectors of rank for each matrices A[i,,]
  d <- dim(A)
  Y <- array(0,d)
  for(i in 1:d[1]){
    Y[i,,] <- projection_rankm(A[i,,], K[i])
  }
  return(Y)
  
}


prod_1 <- function(X, A){
  
  d <- dim(X)
  m <- dim(A)[1]
  Y <- array(1,c(m, d[2], d[3]))
  X_t <- matrix(0,d[1], d[2]*d[3])
  # for(i in 1:d[1]){
  #   X_t[i,]<-c(X[i,,])
  # }
  X_t <- matrix(X, d[1])
  Y_t <- A %*% X_t
  
  # for(i in 1:m){
  #   Y[i,,] <- array(c(Y_t[i,]),d[-1])
  # }
  Y <- array(Y_t,c(m,d[2],d[3]))#NOT SURE OF THAT
  return(Y)
  
}

prod_2_3 <- function(A, B){
  # d1 <- dim(A)
  # d2 <- dim(Q)
  # C <- matrix(0,d1[1], d2[1])
  # for(i in 1:d1[1]){
  #   for(j in 1:d2[1]){
  #     C[i,j] <- crossprod(c(Q[j,,]), c(A[i,,]))
  #   }
  # }
  # A <- matrix(A, nrow = dim(A)[1])
  # Q <- matrix(Q, nrow = dim(Q)[1])
  # C <- tcrossprod(Q,A)
  
  A <- matrix(A, nrow = dim(A)[1])
  B <- matrix(B, nrow = dim(B)[1])
  C <- tcrossprod(A,B)
 
  return(C)
  
  
}

init_algo3 <- function(A,m){
  #Inputs:
  # A: adjacency tensor of dim l x n x n
  # m: number of class of layer
  d <- dim(A)
  A_t <- matrix(0,d[1], d[2]*d[3])
  for(i in 1:d[1]){
    A_t[i,]<-c(A[i,,])
  }
  #A_t <- matrix(A, d[1])
  S <- svds(A_t, m, nu = m)
  W_0 <- S$u
  me <- kmeans(W_0, m, iter.max = 100)
  cent <- me$cluster
  L <- 1:m
  
  Z <- array(0, c(d[1], m))
  for(i in 1:d[1]){
    Z[i,] <- c(L== cent[i])
  }
  
  return(projection_0m(Z))
}


Algo_3 <- function(A, K, iter_max = 150, e = 0.05){
  #Inputs:
  # A: adjacency tensor of dimension Lxnxn
  # K: vector of length M (=number of class of layers) and K_i = number of communities in class i
  M <- length(K)
  L <- dim(A)[1]
  n <- dim(A)[2]
  
  #Initialization:
  W_old <- t(init_algo3(A, M))
  #W_old <- diag(1, M, L)
  #iteration
  i <- 0
  n <- 1
  while((i < iter_max)|(n > e)){
    W <- W_old
    i <- i + 1
    #print(A[1,,]
    Q_1 <- prod_1(A, W_old)
    Q_new <- projection_rankt( Q_1, K)
    #print(Q_new[1,,])
    Z <- prod_2_3(A, Q_new)
    #print(t(Z) %*% Z)
    W_old <- t(projection_0m( prod_2_3(A, Q_new)))
    n <- norm(W_old-W, type = "F")
  }
  
  
  return(list(Q_new, W_old,i))
}


layer_comm_3 <- function(W, m){
  l <- dim(W)[2]
  g <- kmeans(t(W), m)$cluster
  L <- 1:m
  Z <- array(0, c(l, m))
  for(i in 1:l){
    Z[i,] <- c(L== g[i])
  }
  
  return(list(g, Z))
}
nodes_comm_3 <- function(Q, K){
  m <- length(K)
  n <- dim(Q)[2]
  g <- list()
  Z <- list()
  for(i in 1:m){
    Qi <- Q[i,,]
    Ui <- eigen(Qi)$vectors
    Uit <- matrix(0,n, K[i] )
    # for(j in 1:K[i]){
    #   Uit[,j] <- Ui[,j]
    #   
    # }
    Uit <- Ui[,1:K[i]]
    g[[i]] <- kmeans(Uit, K[i])$cluster
    Zi <- array(0, c(n, m))
    L <- 1:m
    for(j in 1:n){
      Zi[j,] <- c(L== g[[i]][j])
    }
    Z[[i]] <- Zi
  }
  return(list(g, Z))
}

