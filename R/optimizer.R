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


parallel_one4_double <- function(...) {
  first <- parallel_one4(...)
  second <- parallel_one4(par=first$par, ...)
}

spg_eps_decreasing <- function(par, control, eps=NULL, ...) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  control$eps <- eps*10
  cat("\033[1;31m")
  zz <- BB::spg(
    par = par,
    ...,
    control = control
  )
  control$eps <- eps
  cat("\033[1;33m")
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  control$eps <- eps/10
  cat("\033[1;32m")
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  
  cat("\033[0m")
  return(zz)
}


spg_eps_decreasing_ftol_dyn <- function(par, control, eps=NULL, noiseseed, fn,
                                        upper,lower,quiet,
                                        ...) {
  cat("starting to compute ftol\n")
  fsd <- sapply(1:30, function(x) fn(theta = par, noiseseed=noiseseed+x, ...)) |> sd()
  control$ftol <-  min(max(fsd/1000,1e-10),1)
  cat("ftol:", control$ftol, "\n")
  if (is.null(eps)) {
    eps <- control$eps
  }
  control$eps <- eps*10
  cat("\033[1;31m")
  zz <- BB::spg(fn = fn,
    par = par,
    ..., noiseseed=noiseseed,
    upper=upper,lower=lower,quiet=quiet,
    control = control
  )
  control$eps <- eps
  cat("\033[1;33m")
  zz <- BB::spg(fn = fn,
    par = zz$par,
    ...,noiseseed=noiseseed,
    upper=upper,lower=lower,quiet=quiet,
    control = control
  )
  control$eps <- eps/10
  cat("\033[1;32m")
  zz <- BB::spg(fn = fn,
    par = zz$par,
    ...,noiseseed=noiseseed,
    upper=upper,lower=lower,quiet=quiet,
    control = control
  )
  
  cat("\033[0m")
  return(zz)
}
spg_eps_decreasing_less <- function(par, control, eps=NULL, ...) {
  if (is.null(eps)) {
    eps <- control$eps
  }
  control$eps <- eps*5
  cat("\033[1;33m")
  zz <- BB::spg(
    par = par,
    ...,
    control = control
  )
  control$eps <- eps/5
  cat("\033[1;32m")
  zz <- BB::spg(
    par = zz$par,
    ...,
    control = control
  )
  
  cat("\033[0m")
  return(zz)
}



parallel_manual_broad_and_fast <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores=25) {
  cat("function: parallel_manual_broad_and_fast")
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
                           cat("inside mcapply\n")
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = 1, noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(150) 
  
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:length(upper)), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 16, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(50) 
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 50, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(10)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.03, triter = 10),
                            prec = precision_factor * 500, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(3) 
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                            prec = precision_factor * 3000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(2)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.005, triter = 10),
                            prec = precision_factor * 8000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(1)
  
  # print(summary(parameters))
  # cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  # cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  # cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
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



parallel_manual_broad_and_fast_mapplymc <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores=25) {
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
      cat(".\n")
      return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(150) 
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
      theta <- c(x1, x2, x3, x4)
      result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                        upper = upper, lower = lower,
                        control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                        prec = precision_factor * 16, noiseseed = noiseseed, ...)
      return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(50) 
  
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
      theta <- c(x1, x2, x3, x4)
      result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                        upper = upper, lower = lower,
                        control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                        prec = precision_factor * 50, noiseseed = noiseseed, ...)
      return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(10)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
      theta <- c(x1, x2, x3, x4)
      result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                        upper = upper, lower = lower,
                        control = list(maximize = FALSE, trace = TRUE, eps = 0.03, triter = 10),
                        prec = precision_factor * 500, noiseseed = noiseseed, ...)
      return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(3) 
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
      theta <- c(x1, x2, x3, x4)
      result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                        upper = upper, lower = lower,
                        control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                        prec = precision_factor * 3000, noiseseed = noiseseed, ...)
      return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(2)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- mcmapply(mc.cores=mc.cores,
    function(x1, x2, x3, x4) {
      theta <- c(x1, x2, x3, x4)
      result <- spg_fun(par = theta, fn = fn, quiet = TRUE,
                        upper = upper, lower = lower,
                        control = list(maximize = FALSE, trace = TRUE, eps = 0.005, triter = 10),
                        prec = precision_factor * 8000, noiseseed = noiseseed, ...)
      return(list(par1 = result$par[1], par2 = result$par[2], par3 = result$par[3], par4 = result$par[4], val = result$value))
    },
    parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |>
    head(1)
  
  # print(summary(parameters))
  # cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  # cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  # cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
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
worked_nicely_once <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores=25) {
  cat("function worked_nicely_once\n")
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
    head(150) 
  
  
  print(summary(parameters))
  
  parameters <- parameters |> head(100)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 50, noiseseed = 1000 * noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(32)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 5),
                            prec = precision_factor * 128, noiseseed = 1000 * noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(8)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.05, triter = 5),
                            prec = precision_factor * 512, noiseseed = 1000 * noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(4)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.02, triter = 5),
                            prec = precision_factor * 2000, noiseseed = 1000 * noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(1)
#parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  # parameters <- parameters |>
  #   rowwise() |>
  #   mutate(
  #     result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
  #                           upper = upper, lower = lower,
  #                           control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 5),
  #                           prec = precision_factor * 8000, noiseseed = 1000 * noiseseed, ...)),
  #     par1 = result$par[1],
  #     par2 = result$par[2],
  #     par3 = result$par[3],
  #     par4 = result$par[4],
  #     val = result$value
  #   ) |>
  #   arrange(val)  |>
  #   head(1)
  # 
  # print(summary(parameters))
  # cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  # cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  # cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
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


parallel_manual_drop_the_last2 <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores=25) {
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
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 50, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(32)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 128, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(8)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.03, triter = 10),
                            prec = precision_factor * 512, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(4)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                            prec = precision_factor * 2000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(1)
  
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  # parameters <- parameters |>
  #   rowwise() |>
  #   mutate(
  #     result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
  #                           upper = upper, lower = lower,
  #                           control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
  #                           prec = precision_factor * 8000, noiseseed = noiseseed, ...)),
  #     par1 = result$par[1],
  #     par2 = result$par[2],
  #     par3 = result$par[3],
  #     par4 = result$par[4],
  #     val = result$value
  #   ) |>
  #select(-result) |>
    #   arrange(val)  |>
  #   head(1)
  # 
  # print(summary(parameters))
  # cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  # cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  # cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
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



parallel_manual_drop_the_last2_flat <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=11,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e5, mc.cores=25) {
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
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 16, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(32)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 128, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(8)
  
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.03, triter = 10),
                            prec = precision_factor * 512, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(4)
  parameters <- rbind(parameters, colmeans(parameters))
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                            prec = precision_factor * 2000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(1)
  
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  # parameters <- parameters |>
  #   rowwise() |>
  #   mutate(
  #     result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
  #                           upper = upper, lower = lower,
  #                           control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
  #                           prec = precision_factor * 8000, noiseseed = noiseseed, ...)),
  #     par1 = result$par[1],
  #     par2 = result$par[2],
  #     par3 = result$par[3],
  #     par4 = result$par[4],
  #     val = result$value
  #   ) |>
#    select(-result) |>
  #   arrange(val)  |>
  #   head(1)
  # 
  # print(summary(parameters))
  # cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  # cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  # cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  # start_time <- Sys.time()
  # 
  # #optimizze again with prec 150000
  # zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
  #             upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
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
parallel_manual <- function(fn, spg_fun=BB::spg, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=15,debug=FALSE,logfn=FALSE, precision_factor=1,   init_cutoff = 1e6, mc.cores=25) {
  tic()
  if (is.null(seed)) {
    noiseseed <- as.integer(runif(1, 1, 1e6))
  } else {
    noiseseed <- seed
  } 
  # draw 1024 points on a grid spanned by lower and upper
  points_per_dimension <- 11
  sequences <- lapply(1:length(lower), function(i) {
    seq(lower[i], upper[i], length.out = points_per_dimension)
  })
  parameters <- expand.grid(sequences)
  
  # loop over grid lines and keep only those with a finite value
  cat("initial grid size:", nrow(parameters), "\n")
  
  cat("1")
  start_time <- Sys.time()
  cat("2")
  colnames(parameters) <- c(paste0("par", 1:length(upper)))

  parameters <- mcmapply(mc.cores=mc.cores,
                         function(x1, x2, x3, x4) {
                           theta <- c(x1, x2, x3, x4)
                           val <- fn(theta, prec = 4, noiseseed = noiseseed, ...)
                           return(list(par1 = x1, par2 = x2, par3 = x3, par4 = x4, val = val))
                         },
                         parameters[,1], parameters[,2], parameters[,3], parameters[,4], SIMPLIFY = F)|> bind_rows()  |> as.data.frame() |> 
    arrange(val)  |> filter(is.finite(val) & val<init_cutoff) |> 
    head(32) 
  
  
  print(summary(parameters))
  
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 1 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  start_time <- Sys.time()
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 16, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(16)
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 2 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  
  # now loop through the 16 points and optimize with spg again
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 128, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(8)
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 3 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.1, triter = 10),
                            prec = precision_factor * 512, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(4)
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 4 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                            prec = precision_factor * 2000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(2) 
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 5 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  parameters <- parameters |>
    rowwise() |>
    mutate(
      result = list(spg_fun(par = c_across(1:4), fn = fn, quiet = TRUE,
                            upper = upper, lower = lower,
                            control = list(maximize = FALSE, trace = TRUE, eps = 0.01, triter = 10),
                            prec = precision_factor * 8000, noiseseed = noiseseed, ...)),
      par1 = result$par[1],
      par2 = result$par[2],
      par3 = result$par[3],
      par4 = result$par[4],
      val = result$value
    ) |>
    select(-result) |>
    arrange(val)  |>
    head(1)
  
  print(summary(parameters))
  cat("best value is:", min(parameters$val), "\nkept", nrow(parameters), "points with finite values, doing next iteration\n")
  cat("best par is:", paste0(parameters[1,1:4] |> unlist(), collapse=", "), "\n")
  cat("round 6 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  start_time <- Sys.time()
  
  #optimizze again with prec 150000
  zz<-spg_fun(par=unlist(parameters[,1:length(lower)]), fn=fn,  quiet=TRUE,
              upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.01, triter=5),
              prec=precision_factor*160000,noiseseed=1000*noiseseed,
              ...)
  
  par<-zz$par
  names(par)<-NULL
  val <- zz$value
  # print par in blue
  cat("\033[34m",par,"\033[0m\n")
  cat("round 7 took ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 0), " minutes\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}

parallel_one4 <- function(fn, lower, upper, seed=NULL, par=NULL, ... ,initialrounds=15,debug=FALSE,logfn=FALSE, precision_factor=1) {
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
                upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.1),
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
              upper=upper,lower=lower,control=list(maximize=FALSE, trace=TRUE, eps=0.001),
              prec=precision_factor*7500,noiseseed=1000*noiseseed, ...
  )
  par<-zz$par
  names(par)<-NULL
  val <- fn(par, ..., prec=precision_factor*10000, noiseseed=noiseseed)
  # print par in blue
  if (debug)
    cat("\033[34m",par,"\033[0m\n")
  
  return(list(par=par,
              val=val,
              tictoc=toc()$callback_msg))   
}


