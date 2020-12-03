# set.seed(1)
# set.seed(round(runif(1)*100000))
# n<-50
# x<-rnorm(n)
# digest::digest(list(x=x))
# y<-0.5*x^2+rnorm(n)*5
# x1<<-NULL
# 
# fit<-lm(y~x)
# print(summary(fit))
# 
# par<-c(coefficients(fit)["(Intercept)"],coefficients(fit)["x"],summary(fit)$sigma  )
# 
# #this is not usefully parallelized, so instead, write a wrapper that does this and
# # takes a set of parameters, a search algorithm and a seed as input
# # runs the command for 1 round and returns and stores the result
# # this needs to be loopable
# 
# #example 
# draw_for_given_par<-function(par,vdata) {
#   BSe<-par[3]*rnorm(n)
#   BSy<-par[1]+par[2]*vdata$x+BSe
#   return(BSy)
# }

bootstrap_parameter_estimate<-function(par,bootstrapseed,optimizer,drawcommand,vdata,...,outcome="y") {
  
  datahash<-digest::digest(vdata)
  optimizerhash<-digest::digest(deparse(optimizer))
  
  dir.create("estimates")
  dir.create(paste0("estimates/",optimizerhash,sep=""))
  dir.create(paste0("estimates/",optimizerhash,"/",datahash,sep=""))

  storedestimates_file<-paste0("estimates/",optimizerhash,"/",datahash,"/boostrapestimates.R",sep="")
  
  callkey<-paste(sep = "",bootstrapseed,digest::digest(c(par,digest::digest(list(...)))))
  
  #open file if exists, use empty list otherwise
  boostrapestimates<-list()
  tryCatch(boostrapestimates<-readRDS(storedestimates_file), error = function(e) {})
  
  if (callkey %in% names(boostrapestimates))
  {
    cat("*")
    return(boostrapestimates[[callkey]])
  }
  set.seed(bootstrapseed)
  #generate fake data
  vdata[[outcome]]<-drawcommand(par,vdata)
  #obtain estimates
  bootest<-optimizer(vdata,...)
  boostrapestimates[[callkey]]<-bootest
  saveRDS(boostrapestimates,storedestimates_file)
  return(bootest)
}
# 
# pars<-NULL
# 
# for (i in 1:1000) {
#   pars<-rbind(pars,c(bootstrap_parameter_estimate(par=par,
#                                                   bootstrapseed=i,
#                                                   drawcommand = draw_for_given_par,
#                                                   vdata=list(x=x),
#                                                   optimizer = function(vdata) {
#                                                     bootfit<-lm(vdata$y ~ vdata$x) 
#                                                     ret<-c(coefficients(bootfit)["(Intercept)"],coefficients(bootfit)["vdata$x"],summary(bootfit)$sigma  )
#                                                     return(ret)
#                                                   }
#   )))
# }
# 
# print(sqrt(diag(var(pars))))
# mean(pars[,2]<0)*2
