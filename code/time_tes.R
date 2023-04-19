library(profvis)
network <- simul_final(100, c(2, 2), c(3,4),2)
Ad <- network[[1]]




profvis({
  r <- Algo_3(network[[1]], c(3, 2))
  
})
