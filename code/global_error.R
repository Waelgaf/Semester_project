


# Error comparing the three approaches (so, one unique class of layer, M=1) ---------------------------

error_123<-function(n, K, L, ro, alpha){
  #Inputs:
  #as usual
  netwe <- alma_generation_M1(n, K, L, 2, ro, alpha )
  A3 <- netwe[[1]]
  A1 <- as.tensor(netwe[[6]])
  L_t <- netwe[[4]]
  true_nod1 <- netwe[[5]][[1]]
  true_nod2 <- netwe[[5]][[2]]
  true_nod <- netwe[[5]]
  #Algorithm 1
  l1 <- which(L_t == 1)
  l2 <- which(L_t == 2)
  #print(l1)
  g11 <- list(Algo_1_new(A3[l1,,], K))
  g12 <- list(Algo_1_new(A3[l2,,], K))
  e1 <- 0.5*(error_nodes_one(true_nod1, g11, K) + error_nodes_one(true_nod2, g12, K)) 
  # print("algo1")
  # print(e1)
  
  #Algorithm 2
  # r <- c(K,K,M)
  # U_init <- InitializationMMSBM(A1, ranks = r )
  # f <- PowerIteration(A1, ranks = r, type="TWIST", U_init, delta1=1000, delta2=1000, max_iter = 25, tol = 1e-05)
  g2 <- nodes_comm_2(A1, L_t, rep(K,2))[[1]]
  # print(g2[[1]])
  e2 <- error_nodes(true_nod, g2, K)
  # print("algo2")
  # print(e2)
  
  #Algorithm 3
  s3 <- Algo_3(A3,rep(K, 2))
  g3 <- nodes_comm_3(s3[[1]], rep(K,2))[[1]]
  e3 <- error_nodes(true_nod, g3, K)
  # print("algo3")
  # print(e3)
  
  return(c(e1,e2,e3))
  
}

# Error comparing the two approaches (two and three) ----------------------


error_12<-function(n, K, L, M, ro, a){
  #Inputs:
  #as usual
  netwe <- alma_generation(n, K, L, M, ro, a )
  A <- netwe[[1]]
  #print(dim(A))
  A1 <- as.tensor(netwe[[6]])
  true_nod <- netwe[[5]]
  true_lay <- netwe[[4]]
  
  #Algorithm 2
  r_1 <- K*M -(M-1) 
  r <- c(r_1,r_1,M)
  U_init <- InitializationMMSBM(A1, ranks = r )
  f <- PowerIteration(A1, ranks = r, type="TWIST", U_init, delta1=1000, delta2=1000, max_iter = 500, tol = 0.01)
  gl2 <- layer_comm_2(f[[2]], M)[[1]]
  gn2 <- nodes_comm_2(A1, gl2, rep(K,M))[[1]]
  en2 <- error_nodes(true_nod, gn2, K)
  el2 <- error_layers(true_lay, gl2, M)
  print(f[[4]]< 500)

  #Algorithm 3
  s3 <- Algo_3(A, rep(K,M), iter_max = 500, e= 0.01)
  gl3 <- layer_comm_3(s3[[2]], M)[[1]]
  gn3 <- nodes_comm_3(s3[[1]], rep(K,M))[[1]]
  en3 <- error_nodes(true_nod, gn3, K)
  el3 <- error_layers(true_lay, gl3, M)
  print(s3[[3]]<500)
  
  return(c(el2, el3, en2, en3))
  
}


