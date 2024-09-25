# Description: This script is used to test the performance of the C++ code to simulate the game 

library(nSGMM)
library(tictoc)
library(igraph)
library(parallel)
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
  simulate_BBP_cpp(n = n,
                   delta0 = -0.2, delta1 = 0.2, sigma = 0.2, 
                   distance = distance, 
                   kinship = kinship, 
                   capacity=capacity,
                   income=income,
                   reps=5, rounds=1000,seed=1)
  
  simulate_BBP_mc(n = n,
                   delta0 = -0.2, delta1 = 0.2, sigma = 0.2, 
                   distance = distance, 
                   kinship = kinship, 
                   capacity=capacity,
                   income=income,
                   reps=5, rounds=1000,seed=1) 
  
  th1<-th1+99.1
  tic()
  for (z in 1:1)
    a<-simulate_BBP_mc(n=n, delta0=th1,delta1=th2,sigma=th3, distance=distance, kinship=kinship, capacity=capacity,income=income, rounds=1000,seed=1,reps=4)  
  toc()
  
  tic()
  for (z in 1:1)
    b<-simulate_BBP_cpp_parallel(n=n, delta0=th1,delta1=th2,sigma=th3, distance=distance, kinship=kinship, capacity=capacity,income=income, rounds=1000,seed=1,reps=4)  
  toc()
  
  tic()
  for (z in 1:1)
    c<-simulate_BBP_cpp(n=n, delta0=th1,delta1=th2,sigma=th3, distance=distance, kinship=kinship, capacity=capacity,income=income, rounds=1000,seed=1,reps=4)  
  toc()
  
  print(a)
  print(b)
  print(c)
  