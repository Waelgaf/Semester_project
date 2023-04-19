library(profvis)
profvis({
  A <- simulation_network(6, 250, 0.5)
  Algo_1(A,3)
  
})

  

