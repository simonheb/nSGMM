library(dplyr)

library(parallel)

simulate_BBP_mc<-function(..., reps, mc.cores=detectCores()-1) {
    finalMatrix <-
    #mclapply(mc.cores=mc.cores,
    lapply(#for debugging
    FUN=function(x) {simulate_BBP_cpp(..., reps=1,indexoffset=x-1)},
    X=1:reps
  ) |> do.call(what=rbind)
  
  if (mean(finalMatrix[,14])<0.20) {
      return(matrix(0,reps,13))
  }
  return(finalMatrix[,1:13]) 
}


simulate_BBP<-function(n,delta0,delta1,sigma,distance,kinship,capacity,income,errors=NULL,seed=1,reps=2,rounds=1000,theta,parallel=FALSE,computeR=FALSE,plotthis=FALSE,kappa.log) {
  oldseed <- .Random.seed
  #cat("repet:",reps,"\n")
  if (any(theta!=c(delta0,delta1,log(sigma),kappatransformation(mean(capacity),log=kappa.log)))) browser() #the duble parameter parameter is not identical
  if (!parallel) {
    finalMatrix<-NULL
    
    for (i in 1:reps) {
      #cat(".",delta0,delta1,sigma,"\n")
      set.seed(seed+i)
      error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
      error <- Rfast::upper_tri.assign(error,rnorm(n*(n-1)/2,sd=sigma))#make symmetric
      error <- Rfast::lower_tri.assign(error,Rfast::lower_tri(t(error)))
      altruism <- 1/(1+exp(-(delta0+delta1*kinship+error)))
      diag(altruism)<-1
      if(mean(Rfast::upper_tri(altruism)>0.999)>0.5) cat("b")
      eq<-equilibrate_and_plot(altruism=altruism,income=income,capacity=capacity,computeR=computeR,computeCPP=!computeR,plotthis = plotthis)
      ret<-compute_moments(1*(eq$transfers>0),kinship,distance,income,theta=c(delta0,delta1,sigma))
      finalMatrix<-rbind(finalMatrix,ret)
    }
  }
  
  .Random.seed<-oldseed
  rnorm(1) #this should update the seed?
  return(finalMatrix)
}
