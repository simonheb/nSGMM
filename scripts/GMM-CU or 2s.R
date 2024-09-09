library(nSGMM)
library(hydroPSO)
setwd("C:/Users/shess/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')

lower<-c(-10,-1,-1,0) 
upper<-c(5,15,3,10) 

############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n <- 25

delta0_DGP<- -1.8
delta1_DGP<- 1
sigma_DGP<-2
kappa_DGP<-2

true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),invkappatransformation(kappa_DGP))

  rseed<-13
  set.seed(rseed)
  # I create the altruism network as a combination of other (random) networks
  kinship <- matrix(rbinom(n*n,1,0.05),nrow=n)
  kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
  diag(kinship)<-1
  
  error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
  error <- Rfast::upper_tri.assign(error,rnorm(n*(n-1)/2,sd=sigma_DGP)) #make symmetric
  error <- Rfast::lower_tri.assign(error,Rfast::lower_tri(t(error)))
  
  altruism <- logistic(delta0_DGP+delta1_DGP*kinship+error)
  diag(altruism)<-1
  
  capacity <- matrix(kappa_DGP,n,n)
  
  income = exp(rnorm(n)+1)

  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity,income=income,plotthis = TRUE)
  
  observed_transfers<-1*(eq$transfers>0)

  
  if (mean(observed_transfers[kinship==1])==1|mean(observed_transfers[kinship==0])==0) 
    cat("kinship perfectly explains failure or success")
  
  # moments to use
  keep <- as.logical(c(0,1,1,1,0,1,0,0,0,0,0,1,0))
  
  vdata<-list(
    kinship=kinship,
    income=income,
    transfers= observed_transfers,
    distance=kinship)
  
  # let's for now assume we know the VCV at the values place
  true_vcv <- moment_distance(theta=true_th,vdata=vdata, keep=keep,prec = 1000)$vcv_full
  
  set.seed(rseed)
  ptm <- Sys.time()
  
  opmizer <- parallel_HydroPSOandSPG
  
  # continuous updating GMM
  result_CU_GMSM <-opmizer(target_function,
                           lower=lower,upper=upper,
                           seed=1, debug=TRUE,
                           vdata=vdata,
                           #vcv=vcv,  
                           keep=keep)
  # 2-step GMM
  result_2S_GMSM_1s <-opmizer(target_function,
                              lower=lower,upper=upper,
                              seed=1, debug=TRUE,
                              vdata=vdata,
                              vcv=diag(length(keep)),  keep=keep)
  fist_vcv <- moment_distance(theta=result_2S_GMSM_1s$par,vdata=vdata, keep=keep,prec = 1000)$vcv_full
  
  result_2S_GMSM <-opmizer(target_function,
                              lower=lower,upper=upper,
                              seed=1, debug=TRUE,
                              vdata=vdata,
                              vcv=fist_vcv,  keep=keep)
  
  # evaluation at true value
  val_CU <- moment_distance(theta=result_CU_GMSM$par,vdata=vdata, keep=keep,prec = 1000,vcv=true_vcv)
  val_2S <- moment_distance(theta=result_2S_GMSM$par,vdata=vdata, keep=keep,prec = 1000,vcv=true_vcv)
  
  