library(nSGMM)
library(tictoc)
library(igraph)
set.seed(12)
  n<-40
  b<-as.matrix(as_adjacency_matrix(erdos.renyi.game(n,0.1+(0.4*runif(1))^2)))
  kinship<-as.matrix(as_adjacency_matrix(erdos.renyi.game(n,(0.5*runif(1))^2)))
  distance<-matrix(rnorm(n*n),nrow=n)
  distance<-distance*t(distance)
  income=rnorm(n)^2
  th1<- -0.2
  th2<- 0.2
  th3<- 0.2
  th4<- 20
  capacity<-matrix(0,n,n)+th4
  zdata<-NULL
  #simulate_BBP_cpp(n, -0.2,0.2,0.2,0.2, distance, kinship, capacity,income,20,1) 
  for(j in 1:100) {
    for(i in seq(1,18,1)) {
    tic()
      for (z in 1:10)
        simulate_BBP_cpp_parallel(n=n, delta0=th1,delta1=th2,sigma=th3, distance=distance, kinship=kinship, capacity=capacity,income=income,theta=c(th1,th2,th3,th4),rounds=1000,seed=1,reps=i)  
    a<-toc()
    zdata<-rbind(zdata,c(i,a$toc-a$tic))
    }
  }
  
  