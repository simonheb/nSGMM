
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(mlrMBO)

setwd("D:/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')

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

lower<-c(-10,-1,-1,-2,-4)
upper<-c(5,15,3,3,2)


print("Somehow i need to think about how the parametrization of the capacity constraint and the parametrization in the 
      estimation can be put together.")
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-40

  delta0_DGP<- -2.4447
  delta1_DGP<- 0.9133
sigma_DGP<-exp(0.4055)
capacity0_DGP<-1.2#65
capacity1_DGP<- -1.3#-0.3

true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity0_DGP,capacity1_DGP)

foo<-NULL
results<-NULL
reps<-1000
for (i in 12:20) {
  source('R/gridsearch.R')
  repseed<-i
  
  cat(i,"=========================================================")
  cat(i,"=========================================================\n")
  rseed<-i
  #set.seed(rseed)
  # #I create the altruism network as a combination of other (random) networks
  lon<-rnorm(n)
  lat<-rnorm(n)
  distance<-as.matrix(dist(cbind(lon,lat),"euclidean",TRUE,TRUE)) #distances between coordinates drawn from a bivariate normal
  distance<-distance/max(distance)*2 #distance in hundresd of meters
  diag(distance)<-0
  
  
  kinship <- matrix(rbinom(n*n,1,0.5),nrow=n)
  kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
  diag(kinship)<-1
  
  error <- matrix(rnorm(n*2,sd=sigma_DGP),nrow=n,ncol=n) 
  error <- upper_tri.assign(error,rnorm(n*(n-1)/2))#make symmetric
  error <- lower_tri.assign(error,lower_tri(t(error)))
  altruism <- 1/(1+exp(-(delta0_DGP+delta1_DGP*kinship+error)))
  capacity <- pmax(capacity0_DGP+capacity1_DGP*distance,0)
  
  diag(altruism)<-1
  
  income = exp(rnorm(n)+1)


  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity,income=income,modmode=21,plotthis = TRUE)
  

  

  
  
  observed_transfers<-1*(eq$transfers>0)
  if (mean(observed_transfers[kinship==1])==1|mean(observed_transfers[kinship==0])==0) cat("XXXXXXXXXXXXXXx")
  
    
  plot(graph_from_adjacency_matrix(observed_transfers))
 
  
  keep<-as.logical(c(1,1,1,0,0,0,1,0,0,1,0,0,1))
  keep_with_degreecorr_geo <-as.logical(c(1,1,1,0,0,0,1,0,0,1,1,0,1))
  keep_with_degreecorr_geo2<-as.logical(c(1,1,1,0,0,0,1,0,1,1,0,0,1))
  
  
  
  
  results<-rbind(results,c(run_1000_new,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  set.seed(repseed)
  ptm <- Sys.time()
  run_1000_new2<-HydroPSOandSPG(g, repfactor=0.25,initialrounds=15, 
                                     lower=lower,upper=upper,  seed=2,
                                     kinship=kinship,
                                     income=income,
                                     transfers= observed_transfers,
                                     distance=distance,
                                     vcv=var(mcpp),  keep=keep_with_degreecorr_geo)
  
  results<-rbind(results,c(run_1000_new2,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  set.seed(repseed)
  ptm <- Sys.time()
  run_1000_new3<-HydroPSOandSPG(g, repfactor=0.25,initialrounds=15, 
                                lower=lower,upper=upper,  seed=2,
                                kinship=kinship,
                                income=income,
                                transfers= observed_transfers,
                                distance=distance,
                                vcv=var(mcpp),  keep=keep_with_degreecorr_geo2)
  
  
  results<-rbind(results,c(run_1000_new3,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  run_1000_new3<-HydroPSOandSPG_fast(g,
                                lower=lower,upper=upper,  seed=2,
                                kinship=kinship,
                                income=income,
                                transfers= observed_transfers,
                                distance=distance,
                                vcv=var(mcpp),  keep=keep_with_degreecorr_geo2)
  
  
  results<-rbind(results,c(run_1000_new3,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  
  print(results )
  saveRDS(results,"a")
  
}
  
  
colmedian<-function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=median, na.rm=na.rm)


