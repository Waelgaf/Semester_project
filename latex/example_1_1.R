B <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.8, 0.3, 0.2, 0.3, 0.8), 3)
Z <- matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1), nrow = 9)
A <- simul_bm(B, Z, 9)
g <- graph_from_adjacency_matrix(A)
cols <- c("blue","red", "yellow")
l <- c(rep(1:3, each = 3))
plot(g, vertex.label=  NA,
     edge.arrow.size = .2, vertex.color=cols[l], edge.color = "grey")
