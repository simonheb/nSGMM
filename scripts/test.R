rm(list=ls())
library(nSGMM)
library(igraph)
set.seed(12)
  n<-20+round(runif(1)*30)
  b<-as.matrix(as_adjacency_matrix(erdos.renyi.game(n,0.1+(0.4*runif(1))^2)))
  kinship<-as.matrix(as_adjacency_matrix(erdos.renyi.game(n,(0.5*runif(1))^2)))
  distance<-matrix(rnorm(n*n),nrow=n)
  distance<-distance*t(distance)
  capacity<-matrix(0,n,n)+9999.1
  income=rnorm(n)^2

  #simulate_BBP_cpp(n, -0.2,0.2,0.2,0.2, distance, kinship, capacity,income,20,1)  
  for(i in 1:10) {
    a<-simulate_BBP_cpp_parallel(n, -0.2,0.2,0.2,0.2, distance, kinship, capacity,income,1000,1)  
  }
  capacity<-matrix(0,n,n)+0.1
  for(i in 1:10) {
    a<-simulate_BBP_cpp_parallel(n, -0.2,0.2,0.2,0.2, distance, kinship, capacity,income,1000,1)  
  }