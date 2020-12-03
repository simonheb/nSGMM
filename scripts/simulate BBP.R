library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(hydroPSO)
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

lower<-c(-10,-1,-1,-2)#,-4)
upper<-c(5,15,3,3)#,2)

#test strategy:
#draw M networks for the estimated parameter values. for each of these, obtain the estimates and recompute the estimates
#this is not usefully parallelized, so instead, write a wrapper that does this and
  # takes a set of parameters, a search algorithm and a seed as input
  # runs the command for 1 round and returns and stores the result
  # this needs to be loopable
#report the standard deviation of each of these estimates




print("Somehow i need to think about how the parametrization of the capacity constraint and the parametrization in the 
      estimation can be put together.")
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-35

  delta0_DGP<- -2
  delta1_DGP<- 1
sigma_DGP<-exp(0.5)
capacity0_DGP<-1.2#65
capacity1_DGP<- NULL#-0.3

true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity0_DGP,capacity1_DGP)

foo<-NULL
results<-NULL




reps<-1000
for (i in 7:129) {
  source('R/gridsearch.R')
  
  
  cat(i,"=========================================================")
  cat(i,"=========================================================\n")
  rseed<-i
  set.seed(rseed)
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
  capacity <- matrix(kappatransformation(capacity0_DGP),n,n) #pmax(capacity0_DGP+capacity1_DGP*distance,0)
  
  diag(altruism)<-1
  
  income = exp(rnorm(n)+1)


  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity,income=income,modmode=21,plotthis = TRUE)
  
  observed_transfers<-1*(eq$transfers>0)

    
    if (mean(observed_transfers[kinship==1])==1|mean(observed_transfers[kinship==0])==0) cat("XXXXXXXXXXXXXXx")
  
    
  plot(graph_from_adjacency_matrix(observed_transfers))
 
  
  keep<-                     as.logical(c(1,1,1,0,0,0,1,0,0,1,0,0,1))
  keep_with_degreecorr_geo2<-as.logical(c(1,1,1,0,0,0,1,0,1,1,0,0,1))
  
  vdata<-list(
    kinship=kinship,
    income=income,
    transfers= observed_transfers,
    distance=distance)
  
  
  set.seed(rseed)
  ptm <- Sys.time()
  run_1000_new3<-HydroPSOandSPG_fast_dep_quick(g,
                                               lower=lower,upper=upper,  seed=1,
                                               vdata=vdata,
                                               initialrounds = 8, debug=TRUE,
                                               vcv=var(mcpp),  keep=keep)
  results<-rbind(results,c(run_1000_new3$par,run_1000_new3$val,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  
  ptm <- Sys.time()
  run_1000_new3<-HydroPSOandSPG_fast_dep_quick(g,
                                               lower=lower,upper=upper,  seed=1,
                                               vdata=vdata,
                                               initialrounds = 10, maxittwo=1,repfactortwo=4, debug=TRUE,
                                               vcv=var(mcpp),  keep=keep)
  results<-rbind(results,c(run_1000_new3$par,run_1000_new3$val,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  
  
  ptm <- Sys.time()
  run_1000_new3<-HydroPSOandSPG_fast_dep_quick(g,
                                     lower=lower,upper=upper,  seed=1, 
                                     vdata=vdata,
                                     maxittwo=1,repfactortwo=4, debug=TRUE,
                                     vcv=var(mcpp),  keep=keep)
  results<-rbind(results,c(run_1000_new3$par,run_1000_new3$val,as.numeric((Sys.time() - ptm),unit="mins"),rseed))
  
  ptm <- Sys.time()
  

  
  colnames(results)<-c("p1","p2","p3","p4","fit","time","round")
  
  print(results )
  saveRDS(results,"a")
  
}



#compare compare fast against slow

g(results[3,1:4], prec=4000,
  kinship=kinship, noiseseed=3,
  income=income, 
  transfers= observed_transfers,
  distance=distance,
  vcv=var(mcpp),  keep=keep)
  
results<-data.frame(results)
results$d1<-results$p1-delta0_DGP
results$d2<-results$p2-delta1_DGP
results$d3<-results$p3-log(sigma_DGP)
results$d4<-results$p4-capacity0_DGP
results$d5<-results$p5-capacity1_DGP

results$type<-c(rep(1:3,4),1,2)
results$d1<-ave(results$d1,results$type,FUN = function(x) {median(abs(x))})
results$d2<-ave(results$d2,results$type,FUN = function(x) {median(abs(x))})
results$d3<-ave(results$d3,results$type,FUN = function(x) {median(abs(x))})
results$d4<-ave(results$d4,results$type,FUN = function(x) {median(abs(x))})
results$d5<-ave(results$d5,results$type,FUN = function(x) {median(abs(x))})
results$tt<-ave(results$time,results$type,FUN = function(x) {median(x)})

