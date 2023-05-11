# library(igraph)
# library(stats)
# library(Matrix)
# library(mc2d)
# library(RSpectra)
# library(combinat)
# library(ggplot2)
library(profvis)
source("Algo3.R")
source("Generation_network.R")
source("error_measures.R")


#Functions which plots the mse error of a generated network (w\ parameters)

MSE_plot_layers_3 <- function(N, k, L, M ){
  
  K <- rep(k, M)
  e <- N
  j <- 1
  for(i in N){
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
    #print(Z_true)
    Z_est <- nodes_comm(A_est[[1]], K)[[1]]
    e[j] <- 0
    
    for(h in 1:M){
      e[j] <- error_nodes_one(Z_true[[h]], Z_est[[h]], k) + e[j]
    }
    j <- j+1
  }
  return(e)
  
}
stoch_block_error_3 <- function(n, k, L, M){
  K <- rep(k, M)
  s <- simul_final(n, K, L, M)
  A_est <- Algo_3(s[[1]], K)
  B_est <- prod_1(A_est[[1]], t(A_est[[2]]))
  heatmap(B_est[1,,])
  return()
  
}

error_algo3 <- function(n, K, L, M){
  s <- simul_final(n, K, L, M)
  A_est <- Algo_3(s[[1]], K)
  #Layer error related
  L_true <- s[[4]]
  L_est <- layer_comm(A_est[[2]], M)[[1]]
  e_l<-  error_layers(L_true, L_est, M)
  
  #Community error related
  Z_true <- s[[5]]
  #print(Z_true)
  Z_est <- nodes_comm(A_est[[1]], K)[[1]]
  e_c <- 0
  
  for(h in 1:M){
    e_c <- error_nodes_one(Z_true[[h]], Z_est[[h]], K) + e_c
  }
  j <- j+1
  
  return(c(e_l,e_c/M))
  
}

error_algo3_new <- function(n, K, L, M){
  s <- simul_final_new(n, K, L, M)
  A_est <- Algo_3(s[[1]], K)
  #Layer error related
  L_true <- s[[4]]
  L_est <- layer_comm(A_est[[2]], M)[[1]]
  e_l<-  error_layers(L_true, L_est, M)
  
  #Community error related
  Z_true <- s[[5]]
  #print(Z_true)
  Z_est <- nodes_comm(A_est[[1]], K)[[1]]
  e_c <- 0
  
  for(h in 1:M){
    e_c <- error_nodes_one(Z_true[[h]], Z_est[[h]], K) + e_c
  }
  j <- j+1
  
  return(c(e_l,e_c/M))
  
}


N <- c(60, 70, 80)
L <- c(10,10, 10)
k <- 3
M <- 3
# profvis({
#   MSE_plot_communities_3(N, k, L, M)
#   
# })
MSE_plot_communities_3(N, k, L, 1)
#MSE_plot_layers_3(N, k, L, M)
# #stoch_block_error_3(60, k, L, M)
#error_algo3_new(100, 3, 40, 3)


