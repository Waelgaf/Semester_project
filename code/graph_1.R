library(igraph)
library(dplyr)
library(ggplot2)
g <- graph_from_adjacency_matrix(Ad[1,,])
cols <- c("blue","red")
q <- x[[1]][1]
plot(g,vertex.shape="none", vertex.label.cex=0.75, edge.color=grey(0.85), edge.width=1, vertex.label.color=cols[q])


