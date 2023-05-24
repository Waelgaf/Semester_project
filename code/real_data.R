library(data.table)
library(igraph)
library(Rcpp)
Rcpp::cppFunction('double mylast(NumericVector x) { int n = x.size(); return x[n-1]; }')
#readLines("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_multiplex.edges", n=10)
dat <- read.table("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_multiplex.edges", sep="")
layers_name <- read.table("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/airport/EUAirTransportation_layers.txt", sep = "")
head(dat)
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
airport_names <- list()
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

#So, Y is the cleared adjacency tensor

#Algo 3 on Y
res <- Algo_3(Y, c(3,3))
W_3 <- res[[2]]
s_r <- eigen(tcrossprod(W_3))$vectors
first_eigen <- t(s_r[1,])%*%tcrossprod(W_3)%*%s_r[1,]
second_eigen <-  t(s_r[2,])%*%tcrossprod(W_3)%*%s_r[2,]
plot(x=first_eigen)