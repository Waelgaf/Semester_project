library(igraph)
library(stats)
library(Matrix)
library(mc2d)
library(RSpectra)
library(combinat)
library(ggplot2)
source("Algo3.R")
source("Generation_network.R")
source("error_measures.R")


#Functions which plots the mse error of a generated network (w\ parameters)

MSE_plot_layers_3 <- function(N, k, L, M ){
  
  K <- rep(k, M)
  e <- N
  j <- 1
  for(i in N){
    print(i)
    s <- simul_final(i, K, L, M)
    A_est <- Algo_3(s[[1]], K)
    L_est <- layer_comm(A_est[[2]], M)[[1]]
    e[j] <- error_layers(s[[4]], L_est, M)
    j <- j + 1
    
  }
  return(e)
}


MSE_plot_communities_3 <- function(N, k, L, M){
  
  K <- rep(k,M)
  e <- N
  j <- 1
  for(i in N){
    s <- simul_final(i, K, L, M)
    A_est <- Algo_3(s[[1]], K)
    Z_true <- s[[5]]
    Z_est <- nodes_comm(A_est[[1]], K)[[1]]
    e[j] <- 0
    
    for(h in 1:M){
      e[j] <- error_nodes_one(Z_true[[h]], Z_est[[h]], k) + e[j]
    }
    j <- j+1
  }
  return(e)
  
}
stoch_block_error_3 <- function(N, k, L, M){
  K <- rep(k, M)
  e <- N
  j <- 1
  for(i in N){
    s <- simul_final(i, K, L, M)
    A_est <- Algo_3(s[[1]], K)
    B_est <- B_est(A_est[[1]], A_est[[2]])
    B_true <- s[[2]]
    # e = Distance between B and \hat{B}
    
  }
  return(e)
  
}


N <- c(60, 70, 80)
L <- c(2,2)
k <- 2
M <- 2

#MSE_plot_communities_3(N, k, L, M)
MSE_plot_layers_3(N, k, L, M)

