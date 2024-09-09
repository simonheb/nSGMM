source('D:/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM/R/gridsearch.R', echo=TRUE)
dims<<-6
ftick<<-0
testfun <- function(theta,noise) {

  
  ftick<<-ftick+1
  gmin<-1:dims
  z<-c(theta-gmin)
  #print(theta)
  z1<-sum(z*z)
  z2<-abs(sin(z))
  return(sum(0.1*z1+t(as.vector(z2))%*%as.vector(z2))+noise*rnorm(1))
}


a<-round(runif(1)*10000)
set.seed(2)
upper<-rep(dims*5,dims)
lower<-rep(-dims*5,dims)
#random_walk_cooling(testfun,theta = c(0,0,0,0) , upper=rep(10,4), minstep = 0.01, lower=rep(-10,4), rounds = 1000, true_theta = 1:dims,
#                    stepsize_rekindling = 10,stepsize_blowup = 2,stepsize_decay = 5)
print(a)
library(hydroPSO)

bounds <- list(lower=lower,upper=upper)

hydroPSO(fn=testfun, lower=bounds$lower, noise=0.00001,
          upper=bounds$upper,
          control=list(maxiter=10000,
            out.with.pbest=TRUE), )
vfdvfdvf


# 
 pkgnames <- c("GA",
               "rgenoud",
               "DEoptim",
               "soma",
               "cmaes",
               "GenSA",
             "pso",
               "NMOF",
               "nloptr",
               "lme4", 
               "hydroPSO",
               "Rmalschains")
 
# ## uncomment to install packages 
# install.packages(pkgnames)
# 
# ## hydromad is not on CRAN
# install.packages(c("zoo", "latticeExtra", "polynom", "car", "Hmisc","reshape"))
# #Install hydromad
# install.packages("hydromad", repos="http://hydromad.catchment.org")

require(hydromad)

for(i in pkgnames) 
  require(i, character.only=TRUE)

require(globalOptTests)


funs <- c("ga", 
#          "genoud",
          "DEoptim",
#          "soma",
#          "cma_es",
          "GenSA", 
#          "psoptim", 
#          "nloptr_crs", 
#          "nloptr_stogo", 
#          "nloptr_d", 
#          "nloptr_d_l", 
#          "nloptr_i", 
#         "optim", 
#          "DEopt", 
          "malschains", 
          "hydroPSO")
#          "SCEoptim",
#          "PSopt"ticks


ans <- array(dim=c( length(funs), 5),
             dimnames=list(funs, NULL))


ticks<-array(0,dim=c( length(funs),2),dimnames = list(funs,NULL))


reps<-10000
  for(j in 1:5) {
    bMat <- cbind(bounds$lower, bounds$upper)
    par <- vector()
    checkDim <- FALSE
    for(ii in 1:length(bounds$lower))
      par <- append(par,runif(1,bounds$lower[ii],bounds$upper[ii])) 
    ## rgenoud 
    
    # out <- try(genoud(fn=testfun, nvars=length(bounds$lower),
    #                   max.generations=10, hard.generation.limit=TRUE,
    #                   pop.size=round(reps/50),
    #                   Domains = matrix(c(bounds$lower,bounds$upper),ncol=2),
    #                   boundary.enforcement=2))
    # if(class(out) == "try-error")
    #   ans["genoud",j] <- NA
    # else
    #   ans["genoud",j] <- max(abs(out$par-1:dims))
    # ticks["genoud",1]<-ftick+ticks["genoud",1]; ftick<<-0
    # 
    # 
    ## ga
    
    minfun<-function(...) {-testfun(...)}
    ftick<<-0
    out <- try(ga(type="real-valued", fitness=minfun, maxiter=1/40*reps,
                  min=bounds$lower, max=bounds$upper, popSize=100))
    if(class(out) == "try-error")
      ans["ga",j] <- NA
    else
      ans["ga",j] <- max(abs(summary(out)$solution-1:dims))
    ticks["ga",1]<-ftick+ticks["ga",1]; ftick<<-0
    
    
    ## DEoptim
    
    ff <- 10*length(par)
    mi <- round(reps/ff)
    out <- try(DEoptim(fn=testfun,lower=bounds$lower, upper=bounds$upper,
                       control=list(itermax=mi)))
    if(class(out) == "try-error")
      ans["DEoptim",j] <- NA
    else
      ans["DEoptim",j] <- max(abs(out$optim$bestmem-1:dims))
    ticks["DEoptim",1]<-ftick+ticks["DEoptim",1]; ftick<<-0
    ## soma
    # 
    # out <- try(soma(costFunction=testfun,bounds=list(min=bounds$lower,
    #                                                 max=bounds$upper),
    #                 options = list(populationSize=round(15*sqrt(reps/5000)),nMigrations=round(30*reps/5000))))
    # if(class(out) == "try-error")
    #   ans["soma",j] <- NA
    # else
    #   ans["soma",j] <- max(abs(out$population[,out$leader]-1:dims))
    # ticks["soma",1]<-ftick+ticks["soma",1]; ftick<<-0
    # ##  cmaes
    # 
    # out <- try(cmaes::cma_es(par=par, fn=testfun,
    #                          lower=bounds$lower, upper=bounds$upper,
    #                          control=list(maxit=1000)))
    # if(class(out) == "try-error")
    #   ans["cma_es",j] <- NA
    # else
    #   ans["cma_es",j] <- out$value
    # 
    #  GenSA
     out <- try(GenSA(fn=testfun,
                      lower=bounds$lower, upper=bounds$upper,
                      control=list(max.call=reps)))
     
     if(class(out) == "try-error")
       ans["GenSA",j] <- NA
     else
       ans["GenSA",j] <- max(abs(out$par-1:dims))
     ticks["GenSA",1]<-ftick+ticks["GenSA",1]; ftick<<-0
    # ##  pso
    
    # out <- try(psoptim(par=par,fn=testfun,
    #                    lower=bounds$lower, upper=bounds$upper,
    #                    control=list(maxit=835)))
    # if(class(out) == "try-error")
    #   ans["psoptim",j] <- NA
    # else
    #   ans["psoptim",j] <- max(abs(out$par-1:dims))
    # ticks["psoptim",1]<-ftick+ticks["psoptim",1]; ftick<<-0
    # 
    ## nloptr crs
    
    #out <- try(crs2lm(x0=par, fn=testfun, 
    #                  lower=bounds$lower, upper=bounds$upper,
    #                  maxeval=10000))
    #if(class(out) == "try-error")
    #  ans["nloptr_crs",j] <- NA
    #else
    #  ans["nloptr_crs",j] <- out$value
    ## nloptr stogo
    
    #out <- try(stogo(x0=par, fn=testfun, 
     #                lower=bounds$lower, upper=bounds$upper,
      #               maxeval=10000))
    #if(class(out) == "try-error")
     # ans["nloptr_stogo",j] <- NA
    #else
     # ans["nloptr_stogo",j] <- out$value
    ## nloptr direct_1
    
    # out <- try(nloptr(x0=par, eval_f=testfun, 
    #                   lb=bounds$lower, ub=bounds$upper,
    #                   opts=list(algorithm="NLOPT_GN_DIRECT_L",
    #                             maxeval=10000)
    #                   ))
    # if(class(out) == "try-error")
    #   ans["nloptr_d_l",j] <- NA
    # else
    #   ans["nloptr_d_l",j] <-  max(abs(out$solution-1:dims))
    # ticks["nloptr_d_l",1]<-ftick+ticks["nloptr_d_l",1]; ftick<<-0
    # ## nloptr direct
    # 
    # out <- try(nloptr(x0=par, eval_f=testfun, 
    #                   lb=bounds$lower, ub=bounds$upper,
    #                   opts=list(algorithm="NLOPT_GN_DIRECT",
    #                             maxeval=10000)
    #                   ))
    # if(class(out) == "try-error")
    #   ans["nloptr_d",j] <- NA
    # else
    #   ans["nloptr_d",j] <-  max(abs(out$solution-1:dims))
    # ticks["nloptr_d",1]<-ftick+ticks["nloptr_d",1]; ftick<<-0
    
    # ## nloptr isres
    # 
    # out <- try(nloptr(x0=par, eval_f=testfun, 
    #                   lb=bounds$lower, ub=bounds$upper,
    #                   opts=list(algorithm="NLOPT_GN_ISRES",
    #                             maxeval=10000)
    #                   ))
    # if(class(out) == "try-error")
    #   ans["nloptr_i",j] <- NA
    # else
    #   ans["nloptr_i",j] <- out$objective 
    ## optim
    
    #out <- try(optim(par=par, fn=testfun, method="SANN", 
    #                 control=list(maxit=10000), 
    #                 ))
    #if(class(out) == "try-error")
    #  ans["optim",j] <- NA
    #else
    #  ans["optim",j] <- out$value
    
    # ## SCEoptim
    # 
    # out <- try(SCEoptim(par=par, FUN=testfun, 
    #                     control=list(maxeval=10000)
    #                     ))
    # if(class(out) == "try-error")
    #   ans["SCEoptim",j] <- NA
    # else
    #   ans["SCEoptim",j] <- out$value
    # 
    ## DEopt
    # 
    # out <- try(DEopt(OF=testfun, algo=list(min=bounds$lower, max=bounds$upper,
    #                                       nG=200/10000*reps, minmaxConstr=TRUE)
    #                  ))
    # if(class(out) == "try-error")
    #   ans["DEopt",j] <- NA
    # else
    #   ans["DEopt",j] <-  max(abs(out$xbest-1:dims))
    # ticks["DEopt",1]<-ftick+ticks["DEopt",1]; ftick<<-0
    # 
    ## PSopt
    
    # out <- try(PSopt(OF=testfun, algo=list(min=bounds$lower, max=bounds$upper,
    #                                       printBar=F, printDetail=F, nG=100)
    #                  ))
    # if(class(out) == "try-error")
    #   ans["PSopt",j] <- NA
    # else
    #   ans["PSopt",j] <- max(abs(out$xbest-1:dims))
    # ticks["PSopt",1]<-ftick+ticks["PSopt",1]; ftick<<-0
    # 
    # 
    ## malschains
    
    
    #goTest1 <- function(par,fnName="testfun",checkDim=F)
     # goTest(par,fnName, checkDim)

     out <- try(malschains(fn=testfun, lower=bounds$lower,
                          upper=bounds$upper,
                          maxEvals=reps))
    

    
    
    if(class(out) == "try-error")
      ans["malschains",j] <- NA
    else
      ans["malschains",j] <- max(abs(out$sol-1:dims))
    ticks["malschains",1]<-ftick+ticks["malschains",1]; ftick<<-0
    ## hydroPSO
    browser()
    out <- try(hydroPSO(fn=testfun, lower=bounds$lower,
                        upper=bounds$upper,
                        control=list(maxfn=reps,out.with.pbest=TRUE), 
                        ))
    if(class(out) == "try-error")
      ans["hydroPSO",j] <- NA
    else
      ans["hydroPSO",j] <-  max(abs(out$par-1:dims))
    ticks["hydroPSO",1]<-ftick+ticks["hydroPSO",1]; ftick<<-0
    
    ## 
    save(ans, file="Results.RData")
    cat("Finished run", j, "under objective function", "testfun", "\n")
  }

print(round(ans,2))
print(round(cbind(ticks,rowMeans(ans)),5))
