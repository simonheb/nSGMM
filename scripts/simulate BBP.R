rm(list=ls())
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
#setwd("D:/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
#source('R/BBP_functions_GMM.R')

 colsums<-colSums
 colmeans<-colMeans
 rowsums<-rowSums
 rowmeans<-rowMeans
 upper_tri<-function(x) {x[upper.tri(x)]}
 lower_tri<-function(x) {x[lower.tri(x)]}
 colMads<-function(x) {colmeans(abs(sweep(x,2,colmeans(x))))}
 rowMads<-function(x) {rowmeans(abs(sweep(x,1,rowmeans(x))))}
 upper_tri.assign<-function(x,y) {x[upper.tri(x)]<-y;return(x)}
 lower_tri.assign<-function(x,y) {x[lower.tri(x)]<-y;return(x)}

 
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

seed<-1#round(runif(1)*100)
set.seed(seed)
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-30
delta0_DGP<- .2
delta1_DGP<- -.2
delta2_DGP<- -.2
sigma_DGP<-.1
capacity_DGP<-1
# #I create the altruism network as a combination of other (random) networks
lon<-rnorm(n)
lat<-rnorm(n)
distance<-as.matrix(dist(cbind(lon,lat),"euclidean",TRUE,TRUE)) #distances between coordinates drawn from a bivariate normal
distance<-distance/max(distance)
diag(distance)<-0

kinship <- matrix(rbinom(n*n,1,0.5),nrow=n)
kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
diag(kinship)<-1

  error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
  error <- upper_tri.assign(error,rnorm(n*(n-1)/2,sd=sigma_DGP))#make symmetric
  error <- lower_tri.assign(error,lower_tri(t(error)))
  altruism <- 1/(1+exp(-(delta0_DGP+delta1_DGP*kinship+delta2_DGP*distance+error)))
  
  diag(altruism)<-1
  #initil incomes are random uniform
    
  income = runif(n)*10
  for(i in 1:20) {
    set.seed(i)
    a<-simulate_BBP_cpp(nrow(kinship),delta0_DGP,delta1_DGP,delta2_DGP,sigma_DGP,distance,kinship,matrix(1,nrow(altruism),ncol(altruism))*capacity_DGP,income,1000)
    b<-simulate_BBP(n=nrow(kinship),delta0=delta0_DGP,delta1=delta1_DGP,delta2=delta2_DGP,sigma=sigma_DGP,distance,kinship,capacity=capacity_DGP,income,reps=1000,parallel=FALSE)
    browser()
  }
#  cdcdcdcdc
 # par(mfrow=c(1,1))
  
  #eq<-equilibrate_and_plot(altruism=altruism,capacity=9,income=income,modmode=21,plotthis = TRUE)

  
  

# 
# #a<-simulate_BBP(nrow(kinship),0.2,-2,-2,2,distance,kinship,99,income,reps=200,modmode=21)

# 
# 
# tictoc::tic()
# cc<-optim(par=c(0,0,0,0.2),upper=rep(1,4),lower=c(-1,-1,-1,0.01),fn=g,transfers=observed_transfers,kinship=kinship,distance=distance,income=income, method = "L-BFGS-B", prec=8, control=list(trace=2,ndeps=c(rep(0.01,4))))
# tictoc::toc()
# tictoc::tic()
# ca<-optim(par=c(0,0,0,0.2),upper=rep(1,4),lower=c(-1,-1,-1,0.01),fn=g,transfers=observed_transfers,kinship=kinship,distance=distance,income=income, method = "L-BFGS-B", prec=8, control=list(trace=2,ndeps=c(rep(0.1,4))))
# tictoc::toc()
# tictoc::tic()
# a<-zoomingGridSearch(g,transfers=observed_transfers,kinship=kinship,distance=distance,income=income, lower=c(-1,-1,-1,0.01),upper=rep(1,4),
#                      stepdepth=7,stepsize=3,stepoverlap=8,stepexpand=1.3,plotit=TRUE,prunepoints = TRUE,pruneratio=0.75)
# tictoc::toc()
# 
# print(g(c(delta0_DGP,delta1_DGP,delta2_DGP,sigma_DGP),transfers=observed_transfers,kinship=kinship,distance=distance,income=income,prec=6))
# print(g(cc$par,transfers=observed_transfers,kinship=kinship,distance=distance,income=income,prec=10))
# print(g(ca$par,transfers=observed_transfers,kinship=kinship,distance=distance,income=income,prec=10))
# print(g(a,transfers=observed_transfers,kinship=kinship,distance=distance,income=income,prec=10))