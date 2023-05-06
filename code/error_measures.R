library(igraph)
library(stats)
library(Matrix)
library(mc2d)
library(combinat)


#Error for one given permutation
error_perm <- function(v, S_true, S_est){
  #Inputs:
  # v: a given permutation
  # S_true: ...
  # S_est: ...
  S_v <- sapply(S_est, function(x) v[x])
  e <- as.vector(S_v) - S_true
  k <- sapply(e, function(x) (x!=0)*1)
  
  return(sum(k))
}


#Missclassification error rate for the layers
error_layers <- function(S_true, S_est, m){
  #Inputs:
  # S_true: true partition of the layer (vector)
  # S_est:estimated partition of the layer (vector)
  L <- length(S_true)
  # L: number of layers
  f <- permn(m, fun = function(x) {error_perm(x,S_true, S_est)})
  print(length(f))
  print(L)
  return(1/L*min(unlist(f)))

}

#Missclassification error rate for a network of one layer for the nodes
error_nodes_one <- function(G_true, G_est, km){
  #Inputs:
  # G_true: true member ship vector of the nodes
  # G_est: estimated membership vector of the nodes
  n <- length(G_true)
  f <- permn(km, fun = function(x){error_perm(x, G_true, G_est)})
  
  return(1/n*min(unlist(f)))
}

