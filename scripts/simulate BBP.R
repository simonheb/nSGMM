rm(list=ls())
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(mlrMBO)

setwd("D:/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')
source('R/gridsearch.R')

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


seed<-2#round(runif(1)*100)
set.seed(seed)
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-30
a<-NULL
delta0_DGP<- -3.1
delta1_DGP<- -1.7
delta2_DGP<- 3.5
sigma_DGP<-0.5
capacity_DGP<-1.4



for (i in 1:15) {
  cat(i,"=========================================================")
  cat(i,"=========================================================")
  cat(i,"=========================================================")
  set.seed(100+i)
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
  
  income = runif(n)^2*10
  
  
  #########emprical vcv####
  simx<-simulate_BBP_cpp_parallel(nrow(kinship),delta0_DGP,delta1_DGP,delta2_DGP,sigma_DGP,
                                  distance,kinship,matrix(capacity_DGP,nrow(kinship),nrow(kinship)),income,1000,1,1000)
  keep<-as.logical(c(1,1,1,1,0,1,1,1))
  
  diff<-sweep(simx,2,colmeans(simx))
  diff<-diff[,keep]
  emprical_vcv<-var(diff)
  diag(emprical_vcv)[which(diag(emprical_vcv)==0)]<-0.00000001 #replace 0 diagonal elements

  
  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity_DGP,income=income,modmode=21,plotthis = TRUE)
  
  
  Sys.sleep(1)
  
  observed_transfers<-1*(eq$transfers>0)
  
  true_th<-c(delta0_DGP,delta1_DGP,delta2_DGP,log(sigma_DGP),capacity_DGP)
  ptm <- proc.time()
  
  trueg<-g(true_th,kinship=kinship,
    income=income,
    transfers=observed_transfers,
    distance=distance,noiseseed=1,prec=1000,verbose=TRUE,vcv=emprical_vcv)
  
  
  
  
  
  ptm <- proc.time()
  run_1000_01_99<-random_walk_cooling(g,c(0,0,0,-1,1),iter=1000,minstep = 0.05,reheatineval=999999,
                                      precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                      lower=c(-5,-5,-5,-5,0.01),upper=c(5,5,5,1,20),
                                      true_theta=c(delta0_DGP,delta1_DGP,delta2_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                      kinship=kinship,
                                      income=income,
                                      transfers=observed_transfers,
                                      distance=distance,noiseseed=1,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_01_99$theta,run_1000_01_99$value,4,(proc.time() - ptm)[3],trueg))
  
  
  ptm <- proc.time()
  run_1000_02<-random_walk_cooling(g,c(0,0,0,-1,1),iter=1000,minstep = 0.2,
                                  precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                  lower=c(-5,-5,-5,-5,0.01),upper=c(5,5,5,1,20),
                                  true_theta=c(delta0_DGP,delta1_DGP,delta2_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                  kinship=kinship,
                                  income=income,
                                  transfers=observed_transfers,
                                  distance=distance,noiseseed=1,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_02$theta,run_1000_02$value,7,(proc.time() - ptm)[3],trueg))
  
  
  
  
}
colmedian<-function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=median, na.rm=na.rm)

g(bottom99$theta,kinship=kinship,
  income=income,
  transfers=observed_transfers,
  distance=distance,noiseseed=3998,prec=500,verbose=TRUE)
g( c( -0.99, 2.46, -3, 0.89, 19.79 ),kinship=kinship,
   income=income,
   transfers=observed_transfers,
   distance=distance,noiseseed=3998,prec=500,verbose=TRUE)
g( c( -0.99-0.1, 2.46-0.1, -3+0.1, 0.89-0.1, 19.79-0.1 ),kinship=kinship,
   income=income,
   transfers=observed_transfers,
   distance=distance,noiseseed=3998,prec=500,verbose=TRUE)

true_th

plot_partial(g,theta= c( -0.99, 2.46, -3, 0.89, 19.79 ),
   kinship=kinship,
  income=income,
  transfers=observed_transfers,
  distance=distance,noiseseed=1,prec=500,param=5,minoffset=-19.6,maxoffset=0,steps=8)



