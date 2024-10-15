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

if (osVersion |> grepl(pattern="Win") |> any())
  mcmapply <- function(..., mc.cores, mc.preschedule) {
    if (osVersion |> grepl(pattern="Win") |> any()) {
      print("windows so using mapply") 
      return(mapply(...))
    }
    else
      return(mcmapply(..., mc.cores=mc.cores, mc.preschedule=mc.preschedule))
  }

spg_eps_decreasing <- function(par, control, eps=NULL, ..., output_id, spg_eps_factor = 10) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  time1 <- Sys.time()
  control$eps <- eps * spg_eps_factor
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
  control$eps <- eps / spg_eps_factor
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
  
  
  cat("spg_eps_decreasing ",output_id,"\t", "*10:", iter2, "in", round(as.numeric(difftime(time2, time1, units = "mins")),2), "mins\t",
      "1:", iter3, "in", round(as.numeric(difftime(time3, time2, units = "mins")),2), "mins\t",
      "*0.1:", iter4, "in", round(as.numeric(difftime(time4, time3, units = "mins")),2), "mins\t",
      zz$value+reduction4+reduction3+reduction2,"=(",reduction2,reduction3,reduction4,")>", zz$value,
      "\n")
  
  return(zz)
}

spg_eps_decreasing_compact <- function(par, control, eps=control$eps, ..., output_id, steps = c(10,1,0.1)) {

  times <- iters <- reductins<-c()
  for (step in steps) {
    starttime <- Sys.time()
    control$eps <- eps * step
    zz <- BB::spg(
      par = par,
      ...,
      control = control
    )
    par <- zz$par
    iters<-c(iters,zz$iter)
    times<-c(times,as.numeric(difftime(Sys.time(), starttime, units = "mins")))
    reductions<-c(reductions,zz$fn.reduction)
  }
  cat("spg_eps_decreasing ",output_id,"\t", "*", step[1], ":", iters[1], "in", times[1], "mins\t",
      step[2], ":", iters[2], "in", times[2], "mins\t",
      "*", step[3], ":", iters[3], "in", times[3], "mins\t",
      zz$value+sum(reductions),"=(",reductions,")>", zz$value,
      "\n")

  return(zz)
}
# 
# shrinking_adaptive_grid <- function(fn, lower, upper,
#                                     seed=1, mc.cores=120,
#                                     start_time=Sys.time(), regularization_lambda=c(1e-3,1e-5,1e-6,1e-7,rep(1e-14,99)),
#                                     startpoint=NULL, depth=20, prec=4, shrinkrate=0.8, radius=(upper-lower)/2,initialrounds=10,
#                                     ...){
#   
#   cat("depth:", depth, "\n")
#   
#   cat("prec:", prec, "\n")
#   cat("regularization_lambda:", regularization_lambda, "\n")
#   cat("how many?", length(regularization_lambda), "\n")
#   
#   if (is.null(seed)) {
#     noiseseed <- as.integer(runif(1, 1, 1e6))
#   } else {
#     noiseseed <- seed
#   } 
#   # draw 1024 points on a grid spanned by lower and upper
#   if (is.null(startpoint)) {
#     cat("function: shrinking_adaptive_grid\n")
#     
#     cat("performance benchmark:")
#     tic()
#     for (i in 1:10)
#       target_function(c(-1,1,1,1), prec = 100, noiseseed = 1, regularization_lambda=0.00001, vdata=vdata, keep=keepsd, maxrounds = 2000, sim_parallel=1 )
#     toc()$callback_msg |> cat()
#     startpoint <- (lower + upper) / 2
#   }
#   
#   grid_lower <- pmax(lower, startpoint - radius)
#   grid_upper <- pmin(upper, startpoint + radius)
#   cat("searching:\n")
#   print(rbind(grid_lower, grid_upper))
#   
#   sequences <- lapply(1:length(lower), function(i) {
#     seq(grid_lower[i], grid_upper[i], length.out = initialrounds)
#   })
#   parameters <- expand.grid(sequences)
#   
#   # loop over grid lines and 
#   
#   colnames(parameters) <- c(paste0("par", 1:length(upper)))
#   
#   if (osVersion |> grepl(pattern="Windows")) {
#     mc.cores <- 1
#     mc.preschedule <- FALSE
#   }
#   parameters <- mclapply(mc.cores=mc.cores, ..., mc.preschedule=FALSE,
#                          FUN = function(theta, ...) {
#                            val <- fn(as.numeric(theta), prec = prec, regularization_lambda = regularization_lambda[1], noiseseed = noiseseed, sim_parallel=FALSE,  ...)
#                            return(c(theta, val = val))
#                          },
#                          X = split(parameters[,1:4],1:nrow(parameters))
#   ) 
#   parameters <- parameters |> bind_rows()  |> as.data.frame() 
#   print(head(parameters,2))
#   parameters <- parameters |> arrange(val)  
#   parameters <- parameters |> head(1)
#   if (depth < 1) {
#     cat("final:\n")
#     print(parameters)
#     cat("took overall:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")),2), "mins\n")
#     return(
#       list(
#         par = parameters[,1:4],
#         value = parameters$val)
#     )
#   }
#   cat("value:\n")
#   print(parameters)
#   shrinking_adaptive_grid(fn = fn, lower=lower, upper=upper, startpoint = parameters[,1:4], start_time=start_time, regularization_lambda=regularization_lambda[2:length(regularization_lambda)], depth = depth-1, prec = pmin(16000,prec*1.5), shrinkrate = shrinkrate, radius = radius*shrinkrate, initialrounds = initialrounds, ...)
#   
# }
# 
# spg_plain <- function(par, control, eps=NULL, ...,output_id) {
#   if (is.null(eps)) {
#     eps <- control$eps
#   }
#   time1 <- Sys.time()
#   control$eps <- eps
#   zz <- BB::spg(
#     par = par,
#     ...,
#     control = control
#   )
#   iter2 <- zz$iter
#   time2 <- Sys.time()
#   zz$step_iter <- c(iter2)
#   zz$step_minutes <- c(as.numeric(difftime(time2, time1, units = "mins")))
#   
#   cat("spg_plain ",output_id,"\t", iter2, "in", round(as.numeric(difftime(time2, time1, units = "mins")),2), "mins\t",
#       zz$value+zz$fn.reduction, "=>", zz$value, "\n")
#   
#   return(zz)
# }

# 
# spg_eps_decreasing_less <- function(par, control, eps=NULL, ...) {
#   if (is.null(eps)) {
#     eps <- control$eps
#   }
#   control$eps <- eps*5
#   time1 <- Sys.time()
#   zz <- BB::spg(
#     par = par,
#     ...,
#     control = control
#   )
#   iter2 <- zz$iter
#   time2 <- Sys.time()
#   control$eps <- eps/5
#   zz <- BB::spg(
#     par = zz$par,
#     ...,
#     control = control
#   )
#   iter3 <- zz$iter
#   time3 <- Sys.time()
#   
#   zz$step_iter <- c(iter2, iter3)
#   zz$step_minutes <- c(as.numeric(difftime(time2, time1, units = "mins")), as.numeric(difftime(time3, time2, units = "mins")))
#   cat("spg_eps_decreasing_less:\t", "*5:", iter2, "in", 
#       round(as.numeric(difftime(time2, time1, units = "mins")),digits = 2), "mins\t",
#       "*0.2:", iter3, "in", round(as.numeric(difftime(time3, time2, units = "mins")),2), "mins\n")
#   
#   return(zz)
# }

sumprogress <- function(round, parameters, start_time) {
  cat("round", round, "took ", floor(10*as.numeric(difftime(Sys.time(), start_time, units = "mins")))/10, " mins\t")
  cat("best value is:", min(parameters$val), "\tkept", nrow(parameters), "pars\n")
  cat("best par is:", paste0(round(parameters[1,1:4]*100)/100 |> unlist(), collapse=", "), "\n")
}


parallel_unified <- function(fn, spg_fun=spg_plain, lower, upper, seed=NULL, par=NULL, ... ,
                             maxit = 1500, 
                             schedule =
                               data.frame(round = c(1,    2,    3,     4,     5,    6,  7),
                                          eps   = c(NA,   0.1,  0.1,  0.1,  0.07,  0.03, 0.01),
                                          keepn = c(240,  50,  25,    12,     6,     3,  1),
                                          precs = c(4,    16,   120, 400,   2000, 4000, 8000),
                                          cutoff_factor = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                                          parallelize_inner = c(F,F,T,T,T,T,T)),
                             regularization = c(1e-3,1e-7,rep(1e-14,nrow(schedule)-2)),
                             initialrounds=14,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5,
                             mc.cores = 120,
                             spg_eps_factor = 10,
                             sim_parallel = 1,
                             mc.preschedule = FALSE
                             ) {
  
  cat("function: parallel_unified\n")

  cat("performance benchmark:")
  tic()
  for (i in 1:10)
    target_function(c(-1,1,1,1), prec = 100, noiseseed = 1, regularization_lambda=0.0001, vdata=vdata, keep=keepsd, maxrounds = 2000, sim_parallel=1 )
  toc()$callback_msg |> cat()
  
  
  start_time_biggi <- Sys.time()
  
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
  
  if (osVersion |> grepl(pattern="Windows")) {
    mc.cores <- 1
    mc.preschedule <- FALSE
  }
  parameters <- mclapply(mc.cores = mc.cores, mc.preschedule = mc.preschedule, ...,
                         FUN = function(theta, ...) {
                           val <- fn(as.numeric(theta), prec = schedule$precs[1], regularization_lambda = regularization[1], noiseseed = noiseseed, ...)
                           return(c(theta, val = val))
                         },
                         X = split(parameters[,1:4],1:nrow(parameters))
                         ) |> 
    bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val)) 
  
  cat("initial grid size after filtering:", nrow(parameters), ", best value is:", min(parameters$val), "\n")
  
  parameters <- parameters |>
    head(schedule$keepn[1]) 
  
  cutoff_val <- min(parameters$val) * schedule$cutoff_factor[1]
  parameters <- parameters |> filter(val < cutoff_val)
  
  
  
  sumprogress(1, parameters, start_time)
  
  
  for (round in 2:nrow(schedule)) {
    start_time <- Sys.time()
    
    if (schedule$parallelize_inner[round]==1) {
      print("parallelizing the inner loop")
      applyfun <- lapply #this means that the mc.cores parameter will be just passed on to FUN by lapply
    } else if (schedule$parallelize_inner[round]==0){
      print("parallelizing the outer loop")
      applyfun <- mclapply #here the mc.cores parameter will not be passed on, instead the the default mc.cores.inner will be used
      mc.cores.inner <- 1
    } else if (schedule$parallelize_inner[round]==2){
      print("parallelizing both ")
      applyfun <- mclapply
      mc.cores.inner <- min(1,floor(mc.cores/nrow(parameters)))
      mc.cores <- nrow(parameters)
    } else {
      stop("parallelize_inner must be 0, 1 or 2")
    }
    parameters <- applyfun(mc.preschedule = mc.preschedule, mc.cores=mc.cores, ...,
                           FUN= function(theta, ..., mc.cores=mc.cores.inner) {
                             result <- spg_fun(par = as.numeric(theta), fn = fn, quiet = TRUE, upper = upper, lower = lower, spg_eps_factor=spg_eps_factor, regularization_lambda = regularization[round], 
                                               control = list(maximize = FALSE, trace = F, eps = schedule$eps[round], triter = 10, maxit = maxit),
                                               prec = precision_factor * schedule$precs[round], noiseseed = noiseseed, ..., output_id = rownames(theta))
                             return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
                           },
                           X = split(parameters[,1:4],1:nrow(parameters))
                           ) |> bind_rows()  |> as.data.frame()
    
    parameters <- parameters |> arrange(val)  |>  head(schedule$keepn[round]) 
    
    cutoff_val <- min(parameters$val) * schedule$cutoff_factor[round]
    parameters <- parameters |> filter(val < cutoff_val)

    
    sumprogress(round, parameters, start_time)
  }
  
  par <- parameters[1,1:4] |> unlist()
  names(par)<-NULL
  val <- min(parameters$val) 
  # print par in blue
  
  cat("\033[34m",par,"\033[0m\n")
  cat("finalvalue:", val, "\n")
  
  cat("took overall:", floor(10*as.numeric(difftime(Sys.time(), start_time_biggi, units = "mins")))/10, " minutes\n")
  return(list(par=par,
              val=val,
              tictoc=floor(10*as.numeric(difftime(Sys.time(), start_time_biggi, units = "mins")))/10))   
}

# 
# parallel_manual_broad_and_fast_mapplymc <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
#                                                     maxit = 1500,
#                                                     initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
#   cat("function: parallel_manual_broad_and_fast_mapplymc")  # print the name of the function that is being exectuted
#   tic()
#   if (is.null(seed)) {
#     noiseseed <- as.integer(runif(1, 1, 1e6))
#   } else {
#     noiseseed <- seed
#   } 
#   # draw 1024 points on a grid spanned by lower and upper
#   sequences <- lapply(1:length(lower), function(i) {
#     seq(lower[i], upper[i], length.out = initialrounds)
#   })
#   parameters <- expand.grid(sequences)
#   
#   # loop over grid lines and keep only those with a finite value
#   cat("initial grid size:", nrow(parameters), "\n")
#   
#   start_time <- Sys.time()
#   
#   colnames(parameters) <- c(paste0("par", 1:length(upper)))
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
#                            return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
#     head(150) 
#   
#   parameters <- rbind(parameters, colmeans(parameters))
#   
#   sumprogress(1,parameters, start_time)
#   
#   start_time <- Sys.time()
#   # now loop through the 16 points and optimize with spg again
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
#                                              upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 16, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |>
#     head(50) 
#   
#   
#   parameters <- rbind(parameters, colmeans(parameters))
#   
#   sumprogress(2,parameters, start_time)
#   
#   
#   start_time <- Sys.time()
#   # now loop through the 16 points and optimize with spg again
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
#                                              upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = FALSE, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 50, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |>
#     head(10)
#   
#   parameters <- rbind(parameters, colmeans(parameters))
#   
#   sumprogress(3,parameters, start_time)
#   
#   start_time <- Sys.time()
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
#                                              upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = FALSE, eps = 0.03, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 500, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |>
#     head(3) 
#   parameters <- rbind(parameters, colmeans(parameters))
#   
#   sumprogress(4,parameters, start_time)
#   
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
#                                              upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = FALSE, eps = 0.01, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 3000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |>
#     head(2)
#   parameters <- rbind(parameters, colmeans(parameters))
#   
#   sumprogress(5, parameters, start_time)
#   
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
#                                              upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = FALSE, eps = 0.005, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 8000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |>
#     head(1)
#   
#   # 
#   #   sumprogress(6,parameters, start_time)
#   # start_time <- Sys.time()
#   # 
#   # #optimizze again with prec 150000
#   # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
#   #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
#   #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
#   #             ...)
#   # 
#   par<-parameters[1,1:4] |> unlist()
#   names(par)<-NULL
#   val <- min(parameters$val) 
#   # print par in blue
#   
#   #sumprogress(7,parameters, start_time)
#   
#   return(list(par=par,
#               val=val,
#               tictoc=toc()$callback_msg))   
# }

# 
# parallel_manual_drop_the_last2 <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
#                                            maxit = 1500,
#                                            initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
#   cat("function: parallel_manual_drop_the_last2")
#   tic()
#   if (is.null(seed)) {
#     noiseseed <- as.integer(runif(1, 1, 1e6))
#   } else {
#     noiseseed <- seed
#   } 
#   # draw 1024 points on a grid spanned by lower and upper
#   sequences <- lapply(1:length(lower), function(i) {
#     seq(lower[i], upper[i], length.out = initialrounds)
#   })
#   parameters <- expand.grid(sequences)
#   
#   # loop over grid lines and keep only those with a finite value
#   cat("initial grid size:", nrow(parameters), "\n")
#   
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
#                            return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> 
#     bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
#     head(100) 
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(1,parameters, start_time)
#   start_time <- Sys.time()
# 
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 50, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val) |> head(32)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(2, parameters, start_time)
#   start_time <- Sys.time()
#   
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 128, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val) |> head(8)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(3, parameters, start_time)
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.03, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 512, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> head(4)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(4, parameters, start_time)
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 2000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> head(2)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(5, parameters, start_time)
#   start_time <- Sys.time()
#   
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 8000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> head(1)
#   sumprogress(6, parameters, start_time) 
#   start_time <- Sys.time()
#   # 
#   # #optimizze again with prec 150000
#   # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
#   #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
#   #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
#   #             ...)
#   # 
#   par<-parameters[1,1:4] |> unlist()
#   names(par)<-NULL
#   val <- min(parameters$val) 
#   # print par in blue
#   cat("\033[34m",par,"\033[0m\n")
#   cat("round 7 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
#   
#   return(list(par=par,
#               val=val,
#               tictoc=toc()$callback_msg))   
# }
# 
# 
# parallel_manual_drop_the_last2_flat <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,
#                                                 maxit = 1500,
#                                                 initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores = 50) {
#   cat("function: parallel_manual_drop_the_last2_flat")
#   tic()
#   if (is.null(seed)) {
#     noiseseed <- as.integer(runif(1, 1, 1e6))
#   } else {
#     noiseseed <- seed
#   } 
#   # draw 1024 points on a grid spanned by lower and upper
#   sequences <- lapply(1:length(lower), function(i) {
#     seq(lower[i], upper[i], length.out = initialrounds)
#   })
#   parameters <- expand.grid(sequences)
#   
#   # loop over grid lines and keep only those with a finite value
#   cat("initial grid size:", nrow(parameters), "\n")
#   
#   start_time <- Sys.time()
#   
#   colnames(parameters) <- c(paste0("par", 1:length(upper)))
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
#                            return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
#     head(100) 
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(1, parameters, start_time)
#   start_time <- Sys.time()
# 
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 16, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val) |> head(32) 
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(2, parameters, start_time)
#   start_time <- Sys.time()
#   
#   
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.1, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 128, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val) |> head(8)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(3, parameters, start_time)
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.03, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 512, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> head(4)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(4, parameters, start_time)
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 2000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val)  |> head(2)
#   parameters <- rbind(parameters, colmeans(parameters))
#   sumprogress(5, parameters, start_time)
#   start_time <- Sys.time()
#   
#   parameters <- mcmapply(mc.cores=mc.cores,
#                          function(x1, x2, x3, x4) {
#                            theta <- c(x1, x2, x3, x4)
#                            result <- spg_fun(par = theta, fn = fn, quiet = TRUE, upper = upper, lower = lower,
#                                              control = list(maximize = FALSE, trace = F, eps = 0.01, triter = 10, maxit = maxit),
#                                              prec = precision_factor * 8000, noiseseed = noiseseed, ...)
#                            return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
#                          },
#                          parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F) |> bind_rows()  |> as.data.frame() |> 
#     arrange(val) |> head(1) 
#   
#   sumprogress(6, parameters, start_time)
#   
#   
#   start_time <- Sys.time()
#   # 
#   # #optimizze again with prec 150000
#   # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
#   #             upper=upper,lower=lower,control=list(maximize=FALSE, trace = FALSE, eps=0.01, triter=5),
#   #             prec=precision_factor*160000,noiseseed=1000*noiseseed,
#   #             ...)
#   # 
#   par<-parameters[1,1:4] |> unlist()
#   names(par)<-NULL
#   val <- min(parameters$val) 
#   # print par in blue
#   cat("\033[34m",par,"\033[0m\n")
#   cat("round 7 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
#   
#   return(list(par=par,
#               val=val,
#               tictoc=toc()$callback_msg))   
# }



parallel_one4 <- function(fn, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=15,debug=FALSE,logfn=FALSE, precision_factor=1) {
  print("parallel_one4")
  tic()
  start_time <- Sys.time()
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
                            prec=precision_factor*7, noiseseed=1000*noiseseed+i, 
                            ...
      )
    else
      z<-hydroPSO::hydroPSO(par=par, fn=lfn, lower=lower, upper=upper,
                            control=list(maxit=100,write2disk=FALSE,verbose=FALSE), 
                            prec=precision_factor*7, noiseseed=1000*noiseseed+i, 
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
  val <- zz$value
  # print par in blue
  cat("final value:", val, "\n")

  cat("\033[34m",par,"\033[0m\n")
  cat("took overall ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}
