library(igraph)
library(stats)
library(Matrix)
library(mc2d)

#We simulate first one single layer from a Stochastic block model
simul_bm <- function(B, Z, n){
  #Input:
  #B: stochastic block matrix
  #Z: membership matrix n x k, w/ k the number of community
  #n :number of individuals
  #Output: Adjacency matrix
  
  
  A <- matrix(0, nrow = n, ncol = n)
  B_u <- Z %*% B %*% t(Z)
  for( i in 1:n){
    for( j in i:n){
      p <- B_u[i,j]
      A[i,j] <- rbern(1,p)
      A[j,i] <- A[i,j]
    }
  }
  return(A-diag(diag(A)))
  
  
}

#We generate here adjacency matrices from a multilayer stochastic block model
simul_mlbm <- function(B, G, n, L){
  #Inputs:
  #B : stochastic block model (list of stochastic block matrix)
  #G : Membership list(for each layer and community): G[[i]] is a n x Ki matrix (~ Z[i,:] in the notes)
  #n : Number of individuals
  #L : vector of size m(=number of layer), membership of class of layer
  #Output: Adjacency tensor
  #L_true: (Layer) Membership vector
  M <- length(B)
  l <- sum(L)
  A <- array(0,c(l,n,n))
  i <- 1
  L_true <- 1:l
  for(m in 1:M){
    for(k in 1:L[m]){
      L_true[i] <- m
      A[i,,] <- simul_bm(B[[m]], G[[m]], n)
      i <- i + 1
    }
  }
  return(list(A, L_true))
}

simul_layer <- function(n, K){
  #Inputs:
  #n : number of individuals
  #K : number of communities
  #Output: Membership matrix
  Z <- matrix(0,n,K)
  Z_true <- 1:n
  for(i in 1:n){
    z <- sample(1:K, 1)
    Z[i,] <- c(z == 1:K)
    Z_true[i] <- z
  }
  return(list(Z,Z_true))
}


simul_final <- function(n, K, L, M = length(L)){
  #Inputs:
  #n: number of individuals
  #K: vector of size M which corresponds to the number of communities in each class of layer
  #L: vectors of layers of size M s.t. li = the number of layer is class i
  #M: number of class of layers
  # Output: 
  # B: Stochastic Blockmodel list
  # A; Adjacency tensor
  # Z: Membership tensor
  # L_true: (Layer) Membership vector
  # G_true:(Communities)list of Membership vector
  l <- sum(L)
  B <- list()
  Z <- list()
  Z_true <- list()
  A <- array(0, c(l,n,n))
  t = 1
  for(i in 1:M){
    sl <- simul_layer(n, K[i])
    Z[[i]] <- sl[[1]]
    Z_true[[i]] <- sl[[2]]
    di <- runif(K[i], 0.9, 1)
    B[[i]] <- diag(di)
    for(j1 in 1:K[i]-1){
      for(j2 in (j1 + 1):K[i]){
        B[[i]][j1,j2] <- runif(1,0,0.5)
        B[[i]][j2,j1] <- B[[i]][j1,j2]
      }
      
    }
    
  }
  S <- simul_mlbm(B, Z, n, L)
  A <- S[[1]]
  L_true <- S[[2]]
  return(list(A, B, Z, L_true, Z_true))
  
}