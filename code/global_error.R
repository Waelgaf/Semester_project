


global_error_2_3 <- function(n, k, L){
  
}

global_error <- function(n, k, L){
  mln <- simul_final_new(n, k, L, 1)
  g_true <- mln[[5]]
  g1_est <- Algo_1(mln[[1]], k)
  g2_est <-
  g3 <- Algo_3(mln[[1]], k)
}