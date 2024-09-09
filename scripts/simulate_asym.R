
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(mlrMBO)

setwd("C:/Users/shess/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
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
for(i in 1:10) {
  simulate_BBP_symmetric_cpp(n,delta0 = delta0_DGP, delta1 = delta1_DGP, sigma = sigma_DGP, distance = distance,kinship = kinship, capacity = matrix(capacity_DGP,n,n),income = income,reps = 400, seed=1,rounds = 333,theta=true_th)
  simulate_BBP_cpp(n,delta0 = delta0_DGP, delta1 = delta1_DGP, sigma = sigma_DGP, distance = distance,kinship = kinship, capacity = matrix(capacity_DGP,n,n),income = income,reps = 400, seed=1,rounds = 333,theta=true_th)
}

seed<-3#round(runif(1)*100)
set.seed(seed)
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-50

  delta0_DGP<- -4.4447
  delta1_DGP<- 1.9133
sigma_DGP<-exp(0.4055)
capacity_DGP<-0.9

true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP)

foo<-NULL

 i<-1

  cat(i,"=========================================================")
  cat(i,"=========================================================\n")
  rseed<-i
  set.seed(rseed)
  # #I create the altruism network as a combination of other (random) networks
  lon<-rnorm(n)
  lat<-rnorm(n)
  distance<-as.matrix(dist(cbind(lon,lat),"euclidean",TRUE,TRUE)) #distances between coordinates drawn from a bivariate normal
  distance<-distance/max(distance)
  diag(distance)<-0
  
  kinship <- matrix(rbinom(n*n,1,0.5),nrow=n)
  kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
  diag(kinship)<-1
  set.seed(as.numeric(Sys.time()))
  error <- matrix(rnorm(n*2,sd=sigma_DGP),nrow=n,ncol=n) 
  altruism <- 1/(1+exp(-(delta0_DGP+delta1_DGP*kinship+error)))
  diag(altruism)<-1
  
  income = exp(rnorm(n)+1)**2

  set.seed(2)
  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity_DGP,income=income,modmode=21,plotthis = TRUE,computeR=F)
  BBP_in_equilibrium_YaT_R(eq$transfers,income = income, altruism = altruism,capacities = matrix(capacity_DGP,n,n),tolerance = 0.001)
  vfvf
  observed_transfers<-1*(eq$transfers>0)
  if (mean(observed_transfers[kinship==1])==1|mean(observed_transfers[kinship==0])==0) cat("XXXXXXXXXXXXXXx")
  
    
  plot(graph_from_adjacency_matrix(observed_transfers))
  #########emprical vcv####
  simx<-simulate_BBP_cpp(nrow(kinship),delta0_DGP,delta1_DGP,sigma_DGP,
                         distance,kinship,matrix(capacity_DGP,nrow(kinship),nrow(kinship)),income,c(delta0_DGP,delta1_DGP,sigma_DGP),1000,1,1000)
  print(cbind(t(compute_moments(observed_transfers,kinship,distance,income,true_th)),compute_moments_cpp(observed_transfers,kinship,distance,income,true_th),colmeans(simx)))

  diff<-sweep(simx,2,colmeans(simx))
  emprical_vcv<-var(diff)
  
  diag(emprical_vcv)[which(diag(emprical_vcv)==0)]<-0.00000001 #replace 0 diagonal elements
  
  
  Sys.sleep(1)
  
  cat("might want to reweight, because this puts an awerfull lot on the moments that identify sigma")
    
  th<-true_th
 

  
  
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, #also use a new moment (c3) to identify the intercept better
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=as.logical(c(1,1,0,1,0,1,1,0,0,0,1,1)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,7201,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,
                  income=income,vcv=emprical_vcv,keep=as.logical(c(1,1,0,1,0,1,1,0,0,0,1,1)),
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  zwei<-run_1000_new
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, #also use a new moment (c3) to identify the intercept better
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=as.logical(c(1,1,1,1,0,0,1,0,0,1,0,0)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,66,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,
                  income=income,vcv=emprical_vcv,keep=as.logical(c(1,1,1,1,0,0,1,0,0,1,0,0)),
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  zwei<-run_1000_new
  
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, #also use a new moment (c3) to identify the intercept better
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=  as.logical(c(1,1,1,0,0,0,1,0,0,1,0,0)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,55,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,
                  income=income,vcv=emprical_vcv,keep=  as.logical(c(1,1,1,0,0,0,1,0,0,1,0,0)),
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  zwei<-run_1000_new
  
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, #how would something neither  path lengths and intermediation work?
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=as.logical(c(1,1,0,1,0,0,1,0,0,0,1,1)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,72011,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,
                  income=income,vcv=emprical_vcv,keep=as.logical(c(1,1,0,1,0,0,1,0,0,0,1,1)),
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  zwei<-run_1000_new
  
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, #use both new moments
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=as.logical(c(1,1,1,1,0,0,1,0,1,0,0,0)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,71,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,  keep=as.logical(c(1,1,1,1,0,0,1,0,1,0,0,0)),
                  income=income,vcv=emprical_vcv,
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  eins<-run_1000_new
  
  
  ptm <- Sys.time()
  run_1000_new<-random_walk_cooling(g,c(0,0,-1,1),iter=1000,minstep = 0.2, 
                                    precschedule=function(iter,maxiter){ceiling(((maxiter-iter)/maxiter)^2*1000+50)},
                                    lower=c(-5,-5,-5,0.01),upper=c(4,4,2,20),
                                    true_theta=c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP),vcv=emprical_vcv,
                                    kinship=kinship, keep=as.logical(c(1,1,0,0,0,0,1,0,0,0,0,1)),
                                    income=income,
                                    transfers=observed_transfers,
                                    distance=distance,momentumdecay = 0.5)
  a<-rbind(a,c(run_1000_new$theta,run_1000_new$value,88,as.numeric((Sys.time() - ptm),unit="mins"),run_1000_new$true_g,rseed,
               g( run_1000_new$theta,kinship=kinship,  keep=as.logical(c(1,1,0,0,0,0,1,0,0,0,0,1)),
                  income=income,vcv=emprical_vcv,
                  transfers=observed_transfers,
                  distance=distance,noiseseed=3998,prec=1000,verbose=TRUE)))
  eins<-run_1000_new
  
  
  
colmedian<-function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=median, na.rm=na.rm)



saveRDS(a,"a")
