library(tictoc)
library(dplyr)
library(parallel)
colMadss<-function(x) {colMeans(abs(sweep(x,2,colMeans(x))))}

colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
plot_partial<-function(theta,param=1,minoffset=-1,maxoffset=1,fun,steps=4,...){
  x<- (0:steps)/steps*(maxoffset-minoffset)+minoffset
  print(x)
  xs<-lapply(x,function(x,theta,par){theta[par] = theta[par] +x;theta},theta=theta,par=param)
  plot(x,lapply(xs,fun,...))
}


spg_eps_decreasing <- function(par, control, eps=NULL, ...) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  time1 <- Sys.time()
  control$eps <- eps*10
  zz <- BB::spg(
    par = par,
    ...,
    control = control
  )
  iter2 <- zz$iter
  time2 <- Sys.time()
  reduction2 <- zz$fn.reduction
  control$eps <- eps
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  iter3 <- zz$iter
  time3 <- Sys.time()
  reduction3 <- zz$fn.reduction
  control$eps <- eps/10
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  iter4 <- zz$iter
  time4 <- Sys.time()
  reduction4 <- zz$fn.reduction
  
  
  zz$step_iter <- c(iter2, iter3,iter4)
  zz$step_minutes <- c(as.numeric(difftime(time2, time1, units = "mins")),
                       as.numeric(difftime(time3, time2, units = "mins")),
                       as.numeric(difftime(time4, time3, units = "mins")))
  
  
  cat("spg_eps_decreasing:\t", "*10:", iter2, "in", round(as.numeric(difftime(time2, time1, units = "mins")),2), "mins\t",
      "1:", iter3, "in", round(as.numeric(difftime(time3, time2, units = "mins")),2), "mins\t",
      "*0.1:", iter4, "in", round(as.numeric(difftime(time4, time3, units = "mins")),2), "mins\t",
      zz$val-reduction4-reduction3,-reduction2,"=(",reduction2,reduction3,reduction4,")>", zz$val,
      "\n")
  
  return(zz)
}



spg_plain <- function(par, control, eps=NULL, ...) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  time1 <- Sys.time()
  control$eps <- eps
  zz <- BB::spg(
    par = par,
    ...,
    control = control
  )
  iter2 <- zz$iter
  time2 <- Sys.time()
  zz$step_iter <- c(iter2)
  zz$step_minutes <- c(as.numeric(difftime(time2, time1, units = "mins")))
  
  cat("spg_plain:\t", iter2, "in", round(as.numeric(difftime(time2, time1, units = "mins")),2), "mins\t",
      zz$val-zz$fn.reduction, "=>", zz$fn.reduction, "\n")
  
  return(zz)
}


spg_eps_decreasing_less <- function(par, control, eps=NULL, ...) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  control$eps <- eps*5
  time1 <- Sys.time()
  zz <- BB::spg(
    par = par,
    ...,
    control = control
  )
  iter2 <- zz$iter
  time2 <- Sys.time()
  control$eps <- eps/5
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  iter3 <- zz$iter
  time3 <- Sys.time()
  
  zz$step_iter <- c(iter2, iter3)
  zz$step_minutes <- c(as.numeric(difftime(time2, time1, units = "mins")), as.numeric(difftime(time3, time2, units = "mins")))
  cat("spg_eps_decreasing_less:\t", "*5:", iter2, "in", 
      round(as.numeric(difftime(time2, time1, units = "mins")),digits = 2), "mins\t",
      "*0.2:", iter3, "in", round(as.numeric(difftime(time3, time2, units = "mins")),2), "mins\n")
  
  return(zz)
}

sumprogress <- function(round, parameters, start_time) {
  cat("round", round, "took ", floor(10*as.numeric(difftime(Sys.time(), start_time, units = "mins")))/10, " minutes\n")
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(round(parameters[1,1:4]*100)/100 |> unlist(), collapse=", "), "\n")
}


parallel_unified <- function(fn, spg_fun=spg_plain, lower, upper, seed=NULL, par=NULL, ... ,
                             maxit = 1500,
                             cutoff_factor = NULL,
                             schedule =
                               data.frame(round = c(1,    2,    3,    4,     5,    6),
                                          eps   = c(NA,  0.1,  0.1, 0.03,  0.01, 0.005),
                                          keepn = c(150, 50,    10,    3,    2,    1),
                                          precs = c(4,    16,   50,  500,  3000, 8000)),
                             initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
  cat("function: parallel_unified\n")
  print(schedule)
  
  tic()
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  } 
  # draw 1024 points on a grid spanned by lower and upper
  sequences <- lapply(1:length(lower), function(i) {
    seq(lower[i], upper[i], length.out = initialrounds)
  })
  parameters <- expand.grid(sequences)
  
  # loop over grid lines and keep only those with a finite value
  cat("initial grid size:", nrow(parameters), "\n")
  
  start_time <- Sys.time()
  
  colnames(parameters) <- c(paste0("par", 1:length(upper)))
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = schedule$precs[1], noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> 
    bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(schedule$keepn[1]) 
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(1, parameters, start_time)
  
  
  for (round in 2:nrow(schedule)) {
    start_time <- Sys.time()
    # now loop through the 16 points and optimize with spg again
    parameters <- mcmapply(mc.cores=mc.cores,
                           function(x1, x2, x3, x4) {
                             theta <- c(x1, x2, x3, x4)
                             result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                               control = list(maximize = FALSE, trace = F, eps = schedule$eps[round], triter = 10, maxit = maxit),
                                               prec = precision_factor * schedule$precs[round], noiseseed = noiseseed, ...)
                             return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                           },
                           parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame()
    
    parameters <- parameters |> arrange(val)  |>  head(schedule$keepn[round]) 
    
    if (!is.null(cutoff_factor)) {
      cutoff_val <- min(parameters$val) * cutoff_factor
      parameters <- parameters |> filter(val < cutoff_val)
    }
    
    if (nrow(parameters)>1) {
      parameters <- rbind(parameters, colmeans(parameters))
    }
    
    sumprogress(round, parameters, start_time)
  }
  
  par <- parameters[1,1:4] |> unlist()
  names(par)<-NULL
  val <- min(parameters$val) 
  # print par in blue
  
  sumprogress("99", parameters, start_time)
  cat("\033[34m",par,"\033[0m\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}


parallel_manual_broad_and_fast_mapplymc <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
                                                    maxit = 1500,
                                                    initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
  cat("function: parallel_manual_broad_and_fast_mapplymc")  # print the name of the function that is being exectuted
  tic()
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  } 
  # draw 1024 points on a grid spanned by lower and upper
  sequences <- lapply(1:length(lower), function(i) {
    seq(lower[i], upper[i], length.out = initialrounds)
  })
  parameters <- expand.grid(sequences)
  
  # loop over grid lines and keep only those with a finite value
  cat("initial grid size:", nrow(parameters), "\n")
  
  start_time <- Sys.time()
  
  colnames(parameters) <- c(paste0("par", 1:length(upper)))
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(150) 
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(1,parameters, start_time)
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                                             upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 16, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(50) 
  
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(2,parameters, start_time)
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                                             upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = FALSE, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 50, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(10)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(3,parameters, start_time)
  
  start_time <- Sys.time()
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                                             upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = FALSE, eps = 0.03, triter = 10, maxit = maxit),
                                             prec = precision_factor * 500, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(3) 
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(4,parameters, start_time)
  
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                                             upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = FALSE, eps = 0.01, triter = 10, maxit = maxit),
                                             prec = precision_factor * 3000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(2)
  parameters <- rbind(parameters, colmeans(parameters))
  
  sumprogress(5, parameters, start_time)
  
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                                             upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = FALSE, eps = 0.005, triter = 10, maxit = maxit),
                                             prec = precision_factor * 8000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(1)
  
  # 
  #   sumprogress(6,parameters, start_time)
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
  #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
  #             ...)
  # 
  par<-parameters[1,1:4] |> unlist()
  names(par)<-NULL
  val <- min(parameters$val) 
  # print par in blue
  
  #sumprogress(7,parameters, start_time)
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}


parallel_manual_drop_the_last2 <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
                                           maxit = 1500,
                                           initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
  cat("function: parallel_manual_drop_the_last2")
  tic()
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  } 
  # draw 1024 points on a grid spanned by lower and upper
  sequences <- lapply(1:length(lower), function(i) {
    seq(lower[i], upper[i], length.out = initialrounds)
  })
  parameters <- expand.grid(sequences)
  
  # loop over grid lines and keep only those with a finite value
  cat("initial grid size:", nrow(parameters), "\n")
  
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> 
    bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(100) 
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(1,parameters, start_time)
  start_time <- Sys.time()

  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 50, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val) |> head(32)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(2, parameters, start_time)
  start_time <- Sys.time()
  
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 128, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val) |> head(8)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(3, parameters, start_time)
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.03, triter = 10, maxit = maxit),
                                             prec = precision_factor * 512, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> head(4)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(4, parameters, start_time)
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
                                             prec = precision_factor * 2000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> head(2)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(5, parameters, start_time)
  start_time <- Sys.time()
  
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
                                             prec = precision_factor * 8000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> head(1)
  sumprogress(6, parameters, start_time) 
  start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
  #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
  #             ...)
  # 
  par<-parameters[1,1:4] |> unlist()
  names(par)<-NULL
  val <- min(parameters$val) 
  # print par in blue
  cat("\033[34m",par,"\033[0m\n")
  cat("round 7 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}


parallel_manual_drop_the_last2_flat <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
                                                maxit = 1500,
                                                initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
  cat("function: parallel_manual_drop_the_last2_flat")
  tic()
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  } 
  # draw 1024 points on a grid spanned by lower and upper
  sequences <- lapply(1:length(lower), function(i) {
    seq(lower[i], upper[i], length.out = initialrounds)
  })
  parameters <- expand.grid(sequences)
  
  # loop over grid lines and keep only those with a finite value
  cat("initial grid size:", nrow(parameters), "\n")
  
  start_time <- Sys.time()
  
  colnames(parameters) <- c(paste0("par", 1:length(upper)))
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(100) 
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(1, parameters, start_time)
  start_time <- Sys.time()

  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 16, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val) |> head(32) 
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(2, parameters, start_time)
  start_time <- Sys.time()
  
  
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
                                             prec = precision_factor * 128, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val) |> head(8)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(3, parameters, start_time)
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.03, triter = 10, maxit = maxit),
                                             prec = precision_factor * 512, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> head(4)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(4, parameters, start_time)
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
                                             prec = precision_factor * 2000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> head(2)
  parameters <- rbind(parameters, colmeans(parameters))
  sumprogress(5, parameters, start_time)
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
                                             control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
                                             prec = precision_factor * 8000, noiseseed = noiseseed, ...)
                           return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
    arrange(val) |> head(1) 
  
  sumprogress(6, parameters, start_time)
  
  
  start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
  #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
  #             ...)
  # 
  par<-parameters[1,1:4] |> unlist()
  names(par)<-NULL
  val <- min(parameters$val) 
  # print par in blue
  cat("\033[34m",par,"\033[0m\n")
  cat("round 7 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}



parallel_one4 <- function(fn, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=15,debug=FALSE,logfn=FALSE, precision_factor=1) {
  print("parallel_one4")
  tic()
  if (logfn) {
    lfn <- function(...) log(fn(...))
  } else {
    lfn <- function(...) fn(...)
  }
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  }
  guess_par <- NULL
  guess_val <- Inf
  for (i in 1:initialrounds) {
    set.seed(noiseseed+i)
    if (is.null(par)) 
      z<-hydroPSO::hydroPSO(fn=lfn, lower=lower, upper=upper,
                            control=list(maxit=100,write2disk=FALSE,verbose=FALSE), 
                            prec=precision_factor*4, noiseseed=1000*noiseseed+i, 
                            ...
      )
    else
      z<-hydroPSO::hydroPSO(par=par, fn=lfn, lower=lower, upper=upper,
                            control=list(maxit=100,write2disk=FALSE,verbose=FALSE), 
                            prec=precision_factor*4, noiseseed=1000*noiseseed+i, 
                            ...
      )
    
    if (debug) {
      #print in color orange, number of it, and the parameters given by hydropsp
      cat("\033[1;33m", "iteration ", i, " z$par: ", z$par, "\n", "\033[0m")
    }
    
    
    zz<-BB::spg(par=z$par, fn=lfn,  quiet=TRUE,
                upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.1),
                prec=precision_factor*300,noiseseed=1000*noiseseed,
                ...
    )
    if (debug) {
      #print in color orange, number of it, and the parameters given by spg
      cat("\033[1;33m", "iteration ", i, " zz$par: ", zz$par, "\n", "\033[0m")
    }
    
    if (zz$value < guess_val) {
      guess_par <- zz$par
      guess_val <- zz$value
    }
  }
  
  print("final_spg:")
  zz<-BB::spg(par=guess_par, fn=lfn, 
              upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.001),
              prec=precision_factor*7500,noiseseed=1000*noiseseed, ...
  )
  par<-zz$par
  names(par)<-NULL
  val <- zz$val
  # print par in blue
  if (debug)
    cat("\033[34m",par,"\033[0m\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}
