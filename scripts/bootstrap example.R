rm(list=ls())
set.seed(1)
set.seed(round(runif(1)*100000))
n<-50
x<-rnorm(n)
digest::digest(list(x=x))
y<-0.5*x^2+rnorm(n)*5
x1<<-NULL

fit<-lm(y~x)
par<-c(coefficients(fit)["(Intercept)"],coefficients(fit)["x"],summary(fit)$sigma  )

bootcoeffs<-NULL
for (i in 1:1000) {
  BSy<-par[1]+par[2]*x+par[3]*rnorm(n)
  bootfit<-lm(BSy ~ x)
  bootcoeffs<-c(bootcoeffs,coefficients(bootfit)["x"])
}
print(round(coefficients(summary(fit)),3))
cat("x            ",round(mean(bootcoeffs),3),"   ",round(sd(bootcoeffs),3),"\n")



#this is not usefully parallelized, so instead, write a wrapper that does this and
# takes a set of parameters, a search algorithm and a seed as input
# runs the command for 1 round and returns and stores the result
# this needs to be loopable

draw_for_given_par<-function(par,vdata) {
  BSe<-par[3]*rnorm(n)
  BSy<-par[1]+par[2]*vdata$x+BSe
  return(BSy)
}

bootstrap_parameter_estimate<-function(par,bootstrapseed,optimizer,drawcommand,vdata,...) {
  boostrapestimates<-list()
  tryCatch(boostrapestimates<-readRDS("boostarapdestimates"), error = function(e) {})
  callkey<-paste(sep = "",bootstrapseed,digest::digest(c(par,digest::digest(vdata),digest::digest(deparse(optimizer)),digest::digest(list(...)))))
  
  if (callkey %in% names(boostrapestimates))
  {
    cat("*")
    return(boostrapestimates[[callkey]])
  }
  set.seed(bootstrapseed)
  #generate fake data
  observed<-drawcommand(par,vdata)
  #obtain estimates
  bootest<-optimizer(booty=observed,vdata,...)
  boostrapestimates[[callkey]]<-bootest
  saveRDS(boostrapestimates,"boostarapdestimates")
  return(bootest)
}

pars<-NULL

for (i in 1:1000) {
  pars<-rbind(pars,c(bootstrap_parameter_estimate(par=par,
                                                  bootstrapseed=i,
                                                  drawcommand = draw_for_given_par,
                                                  vdata=list(x=x),
                                                  optimizer = function(booty,vdata) {
                                                    bootfit<-lm(booty ~ vdata$x) 
                                                    ret<-c(coefficients(bootfit)["(Intercept)"],coefficients(bootfit)["vdata$x"],summary(bootfit)$sigma  )
                                                    return(ret)
                                                  }
  )))
}

print(sqrt(diag(var(pars))))
mean(pars[,2]<0)*2
