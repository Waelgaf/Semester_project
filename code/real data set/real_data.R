set.seed(44)
library(data.table)
library(igraph)
library(Rcpp)
Rcpp::cppFunction('double mylast(NumericVector x) { int n = x.size(); return x[n-1]; }')
dat <- read.table("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_multiplex.edges", sep="")
layers_name <- read.table("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_layers.txt", sep = "")
nodes_name <- read.table("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_nodes.txt", sep = "")
d <- dat$V1
l <- mylast(d)
A <- array(0, c(l, 450,450))
k <- as.matrix(dat)
r <- dim(k)
for(j in 1:r[1]){
  i <- k[j,]
  i1 <- i[2]
  i2 <- i[3]
  m <- i[1]
  A[m,i1, i2] <- 1
  A[m, i2, i1] <- 1
}
g <- list()
i <- 0
deg <- rep(NA, l)
Y <- array(NA, c(12, 450, 450))
companies <- rep(NA,12)
for(j in 1:l){
  g[[j]] <- graph_from_adjacency_matrix(A[j,,], mode= "undirected")
  d <- mean(degree(g[[j]]))
  if(d > 0.35){
    i <- i+1
    companies[i] <- j
    Y[i,,]<- A[j,,]
  }
  
}

companies_names <- layers_name[companies,]


#Remove the (non-visited airports)
k <- 0
clean <- c()
A_f <- array(NA, c(12, 242, 242))
for(i in 1:450){
  if(mean(Y[,i,]) < 0.00045){
    clean <- append(clean, i)
    k <- k+1
  }
}

A_f <- Y[,-clean, -clean]

airports <- nodes_name$V2[-clean]
#So, A_f is the cleared adjacency tensor
#Define the spectral distance

msqrt_d <- function(d){
  d1 <- c()
  for(i in d){
    if(i == 0){
      d1 <- append(d1,0)
    }else{
      d1 <- append(d1, i^(-1/2))
    }
  }
  return(diag(d1))
}

normalized_Lapl <- function(A){
  d <- apply(A, 1, sum)
  D <- diag(d)
  D_12 <- msqrt_d(d)
  return(D_12%*%(D-A)%*%D_12)
}
spectr_dist <- function(A, B){
  L1 <- normalized_Lapl(A)
  L2 <- normalized_Lapl(B)
  l1 <- eigen(L1, symmetric = TRUE, only.values = TRUE)$values
  l2 <- eigen(L2, symmetric = TRUE, only.values = TRUE)$values
  return(sum(sapply(l1-l2, function(x)abs(x))))
}
#Dendrogram of the layers of A_f
A_lf <- list()
for(i in 1:12){
  A_lf <- c(A_lf, list(A_f[i,,]))
}
distanc <- array(NA, c(12,12))
for(i in 1:12){
  for(j in i:12){
    distanc[i,j] <- spectr_dist(A_lf[[i]], A_lf[[j]])
    distanc[j,i] <- distanc[i,j]
  }
}
distanc <- as.data.frame(distanc, row.names = companies_names$V2)
plot(hclust(as.dist(distanc)), hang = -1, cex = 1)

#It makes sense to consider M=3
#Now, we want to choose the number of communities K = (K,K,K)
n <- 242
l <-12
m <-3
  
A_tf <- array(NA, c(n,n,l))
for(i in 1:l){
  print(i)
  A_tf[,,i] <- A_f[i,,]
}
A_tf <- as.tensor(A_tf)
K <- c(2,3,4,5,6)
tot_var2 <- rep(NA, length(K))
tot_var3 <- rep(NA, length(K))
j <- 0
for(k in K){
  q <- k*m -(m-1)
  r <- c(q,q,m)
  j <- j + 1
  print("a")
  #Algorithm 2
  U_init <- InitializationMMSBM(A_tf, ranks = r )
  s2 <- PowerIteration(A_tf, ranks = r, type="TWIST", U_init, delta1=1000, delta2=1000, max_iter = 500, tol = 0.001)
  gl2 <- layer_comm_2(s2[[2]], m)[[1]]
  gn2 <- nodes_comm_2_real_data(A_tf, gl2, rep(k,m))
  tot_var2[j] <- gn2[[3]]
  print("b")
  #Algorithm 3
  s3 <- Algo_3(A_f, rep(k,m), iter_max = 500, e= 0.001)
  gl3 <- layer_comm_3(s3[[2]], m)[[1]]
  gn3 <- nodes_comm_2_real_data(A_tf, gl3, rep(k,m))
  tot_var3[j] <- gn3[[3]]
  print("e")
}
#Elbow method on graph g
g <- ggplot(data = NULL, aes(x=K))+ geom_line(aes(y = tot_var2), colour = "blue", size = 1) + geom_line(aes(y = tot_var3), colour = "red", size = 1)+ylab("sum of squares")
g + theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
#So, k = 4, (m=3)
k <- 4
m <- 3
#Algo 2
q <- k*m -(m-1)
r <- c(q,q,m)
U_init <- InitializationMMSBM(A_tf, ranks = r )
s2 <- PowerIteration(A_tf, ranks = r, type="TWIST", U_init, delta1=1000, delta2=1000, max_iter = 500, tol = 0.001)
gl2 <- layer_comm_2(s2[[2]], m)[[1]]
gn2 <- nodes_comm_2(A_tf, gl2, rep(k,m))[[1]]

glay2 <- ggplot(NULL, aes(x = s2[[2]][,2], y = s2[[2]][,3], label = companies_names$V2)) +geom_text(check_overlap = TRUE, aes(color = factor(gl2)))
glay2 + xlab("second component") + ylab("third component") + theme(legend.position = "none")+theme(axis.text=element_text(size=12),
                                                                                                   axis.title=element_text(size=12))



#Algo3
s3 <- Algo_3(A_f, rep(k,m), iter_max = 500, e= 0.001)
gl3 <- layer_comm_3(s3[[2]], m)[[1]]
gn3 <- nodes_comm_3(s3[[1]], rep(k,m))[[1]]

glay3 <- ggplot(NULL, aes(x = s3[[2]][2,], y = s3[[2]][3,], label = companies_names$V2)) +geom_text(check_overlap = FALSE, aes(color = factor(gl3)))
glay3 + xlab("second component") + ylab("third component") + theme(legend.position = "none")+theme(axis.text=element_text(size=12),
                                                                                                   axis.title=element_text(size=12))
cols <- c(  "seagreen4", "chocolate1", "mediumblue",  "darkviolet")

g2_nodes_ryanair <- gn2[[2]]
g3_nodes_ryanair <- gn3[[2]]
ryanair_graph <- graph_from_adjacency_matrix(A_f[2,,], mode= "undirected")
isolated <- which(degree(ryanair_graph)==0)
#ryanair_graph <- delete.vertices(ryanair_graph, isolated)
plot(ryanair_graph, vertex.size = degree(ryanair_graph)/20, vertex.label=  airports, vertex.label.cex = .8, vertex.label.color = cols[g2_nodes_ryanair],
edge.arrow.size = .02, vertex.color= "snow3", edge.color = "gray79")

B <- as.data.frame(A_f[2,-isolated, -isolated], col.names = airports[-isolated], row.names = airports[-isolated])
heatmap(A_f[2, -isolated, -isolated], Colv = NA, Rowv = NA, symm = TRUE, labRow = airports[-isolated], labCol = airports[-isolated] , col =  heat.colors(5))
