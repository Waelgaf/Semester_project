library(profvis)
#profvis({
  network <- simul_final(250, c(2, 2), c(3,4),2)
  Ad <- network[[1]]
  r <- Algo_3(Ad, c(3, 2))
  y <- layer_comm(r[[2]], 2)
  x <- nodes_comm(r[[1]], c(2,2))
#})
