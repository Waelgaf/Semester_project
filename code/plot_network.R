library(igraph)
library(stats)
library(Matrix)
library(mc2d)

source("Algo1.R")
source("Algo3.R")
source("Generation_network.R")

plot_network <- function(A ,K, k){
  #Inputs:
  # A: Adjacency matrix
  # K: Membership vector
  # k: number of communities
  g <- graph_from_adjacency_matrix(A)
  
  cols <- rainbow(k)
  
  plot(g,vertex.size = degree(g)/5, vertex.label=  NA,
      edge.arrow.size = .2, vertex.color=cols[q], edge.color = "grey")
  
}

plot_networks <- function(A, L, K){
  #Inputs:
  # A: Adjacency tensor
  # L: Membership layer class
  # K: Membership community (list of vector)
  # k: number of communities in each class of layer (vector)
  l <- dim(A)[1]
  
  for(i in 1:l){
    plot_network(A[i,,], unlist(K[[i]]), k[i])
  }
  
  
}