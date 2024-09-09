
rm(list=ls())
library(nSGMM)
library(hydroPSO)
library(tictoc)
library(Rfast)
library(igraph)
setwd("C:/Users/shess/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')
source("R/optimizer.R")




# moments to use

keepsd   <-   as.logical(c(1,1,1,1,1,1,0,0,0,0,0,0,0))
keepsdNib <-  as.logical(c(1,1,0,1,1,1,0,0,0,0,0,0,0))



lower<-c(-11,-11,-5,0) 
upper<-c(11,11,5,10)

vsdata<-readRDS("../gambiafitting/DataforBBP")
transfernet<-"m4m6m7"

############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 

set.seed(1)
# take any village from vsdata
vill <- sample(1:length(vsdata),1)
vdata <- vsdata[[vill]]
n <- length(vdata$income)
print(n) # 52 69 56 34 43

delta0_DGP <- -7
delta1_DGP <- 5
sigma_DGP <- 4
kappa_DGP <- 5


# This gives a nice difficult one
# 
# lower<-c(-10,-10,-5,0) 
# upper<-c(10,10,5,10)
# 
# vsdata<-readRDS("../gambiafitting/DataforBBP")
# transfernet<-"m4m6m7"
# 
# ############# Set Simulation Parameters and Draw Random Altruism Network#############
# #population size 
# 
# set.seed(543532411)
# # take any village from vsdata
# vill <- sample(1:length(vsdata),1)
# vdata <- vsdata[[vill]]
# n <- length(vdata$income)
# print(n) # 32 67 60 37 39 52 69 56 34
# 
# set.seed(13)
# delta0_DGP <- -6
# delta1_DGP <- -3
# sigma_DGP <- 6
# kappa_DGP <- 25


true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),invkappatransformation(kappa_DGP)) #-0.5 1.5 -0.693 7.408

# I create the altruism network as a combination of other (random) networks
kinship <- vdata$m8am8bm8c
#kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
#diag(kinship)<-1




income = vdata$income+1
#set.seed(as.numeric(Sys.time()))
error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
error <- Rfast::upper_tri.assign(error,rnorm(n*(n-1)/2,sd=sigma_DGP)) #make symmetric
error <- Rfast::lower_tri.assign(error,Rfast::lower_tri(t(error)))

altruism <- logistic(delta0_DGP+delta1_DGP*kinship+error)


diag(altruism)<-1


capacity <- matrix(kappa_DGP,n,n) 

hist(altruism[!diag(n)])

eq <- equilibrate_and_plot(altruism=altruism,capacity=capacity,income=income,plotthis = TRUE)


# binarize
observed_transfers<-1*(eq$transfers>0)
 


if (mean(observed_transfers[kinship==1])==1|mean(observed_transfers[kinship==0])==0) 
  cat("kinship perfectly explains failure or success")

vdata<-list(
  kinship=kinship,
  income=income,
  transfers= observed_transfers,
  distance=kinship)


set.seed(12345)
# lets look at different variangs of sdNiB
# c <- 
#   rbind(t(unlist(run_1000_new1_cugmm_sdNiB)),
#         t(unlist(run_1000_new1_cugmm_sdNiB_regul)),
#         t(unlist(run_1000_new1_cugmm_sdNiB_v1)),
#         t(unlist(run_1000_new1_cugmm_sdNiB_regul_v1))
#   ) |> as.data.frame()
# c$time <- gsub("[^0-9.]*", replacement = "",  c$tictoc)
# c$tictoc <- NULL
# c$type <- c(1,2,1,2)
# c <- as.data.frame(lapply(c, as.numeric))
# 
# library(ggplot2)
# ggplot(c,aes(y=(val),x=time,label=type))+
#   geom_point()+
#   geom_text(nudge_y=0.01)+
#   theme_minimal()
# 
# 
# jhkghkhjkhjk
## Check where seeds go wrong:





# nochmal_ein_versuch_full <- parallel_manual(target_function,
#                                        lower=lower,upper=upper,
#                                        seed=1, 
#                                        vdata=vdata,
#                                        regularization_lambda = 1e-14,
#                                        keep=keepsd)


nochmal_ein_versuch_broad_and_fast <- parallel_manual_broad_and_fast(target_function,
                                                                     spg_fun = spg_eps_decreasing,
                                                                     initialrounds = 8,
                                                                     lower=lower,upper=upper,
                                                                     seed=1, 
                                                                     vdata=vdata,
                                                                     regularization_lambda = 1e-14,
                                                                     keep=keepsd)

nochmal_ein_versuch_shortdec8 <- parallel_manual_drop_the_last2(target_function,
                                                                spg_fun = spg_eps_decreasing,
                                                                initialrounds = 8,
                                                                lower=lower,upper=upper,
                                                                seed=1, 
                                                                vdata=vdata,
                                                                regularization_lambda = 1e-14,
                                                                keep=keepsd)






nochmal_ein_versuch_shortdec8_less <- parallel_manual_drop_the_last2(target_function,
                                                                spg_fun = spg_eps_decreasing_less,
                                                                initialrounds = 8,
                                                                lower=lower,upper=upper,
                                                                seed=1, 
                                                                vdata=vdata,
                                                                regularization_lambda = 1e-14,
                                                                keep=keepsd)

nochmal_ein_versuch_shortdec8_ftol_dyn <- parallel_manual_drop_the_last2(target_function,
                                                                spg_fun = spg_eps_decreasing_ftol_dyn,
                                                                initialrounds = 8,
                                                                lower=lower,upper=upper,
                                                                seed=1, 
                                                                vdata=vdata,
                                                                regularization_lambda = 1e-14,
                                                                keep=keepsd)


nochmal_ein_versuch_worked_nicely_once <- worked_nicely_once(target_function,
                                                             spg_fun = spg_eps_decreasing,
                                                             initialrounds = 8,
                                                             lower=lower,upper=upper,
                                                             seed=1, 
                                                             vdata=vdata,
                                                             regularization_lambda = 1e-14,
                                                             keep=keepsd)

  

nochmal_ein_versuch_short8 <- parallel_manual_drop_the_last2(target_function,
                                                             initialrounds = 8,
                                                             lower=lower,upper=upper,
                                                             seed=1, 
                                                             vdata=vdata,
                                                             regularization_lambda = 1e-14,
                                                             keep=keepsd)
nochmal_ein_versuch_short <- parallel_manual_drop_the_last2(target_function,
                                                            lower=lower,upper=upper,
                                                            seed=1, 
                                                            vdata=vdata,
                                                            regularization_lambda = 1e-14,
                                                            keep=keepsd)

nochmal_ein_versuch_shortdec <- parallel_manual_drop_the_last2(target_function,
                                                               spg_fun = spg_eps_decreasing,
                                                               lower=lower,upper=upper,
                                                               seed=1, 
                                                               vdata=vdata,
                                                               regularization_lambda = 1e-14,
                                                               keep=keepsd)
nochmal_ein_versuch_shortdec_less <- parallel_manual_drop_the_last2(target_function,
                                                               spg_fun = spg_eps_decreasing_less,
                                                               lower=lower,upper=upper,
                                                               seed=1, 
                                                               vdata=vdata,
                                                               regularization_lambda = 1e-14,
                                                               keep=keepsd)

#nochmal_ein_versuch_full <- readRDS("nev.rds")

# 
# run_1000_new1_cugmm_sdNiB <- parallel_one4(target_function,
#                                            lower=lower,upper=upper,
#                                            seed=1, 
#                                            vdata=vdata, 
#                                            regularization_lambda = 1e-14,
#                                            keep=keepsdNib)

run_1000_new1_cugmm_sd <- parallel_one4(target_function,
                                        lower=lower,upper=upper,
                                        seed=1, 
                                        vdata=vdata,
                                        regularization_lambda = 1e-14,
                                        keep=keepsd)


nochmal_ein_versuch_shortdec8_flat <- parallel_manual_drop_the_last2_flat(target_function,
                                                                          spg_fun = spg_eps_decreasing,
                                                                          initialrounds = 8,
                                                                          lower=lower,upper=upper,
                                                                          seed=1, 
                                                                          vdata=vdata,
                                                                          regularization_lambda = 1e-14,
                                                                          keep=keepsd)

bgbfgb
#nochmal_ein_versuch <- run_1000_new1_cugmm_sd
# drop anything as soon as I reach 10 strikes
a <- 
  rbind(#t(unlist(run_1000_new1_cugmm_sd)),      # ||||| ||||| |
        #t(unlist(run_1000_new1_cugmm_sdNiB)),   # ||||| ||||| ||||| |||||
        #t(unlist(nochmal_ein_versuch)),         # ||||| |||
        t(unlist(nochmal_ein_versuch_short)),   # ||||| |||
        t(unlist(nochmal_ein_versuch_shortdec))   
  ) |> as.data.frame()
a$time <- gsub("[^0-9.]*", replacement = "",  a$tictoc) 
a$tictoc <- NULL
# make all entries in a numeric
a <- as.data.frame(lapply(a, as.numeric))
a$type <- c(1:nrow(a))
#compute new function value anew for each of the rows
a$val_new <- 0

library(ggplot2)
ggplot(a,aes(y=log(val),x=time,label=type))+
  geom_point()+
  geom_text(nudge_y=0.0015)+
  theme_minimal()


se_theta<-c(3,1.8,0.3,2.3)
b <- rbind(
  true_th,
  a[,1:4]
)

deviations_from_true<-t(t(b)-true_th)
standardized_deviations <- (abs(t(t(deviations_from_true)/ se_theta)))
standardized_deviations <- standardized_deviations[2:nrow(standardized_deviations),] |> as.data.frame()

newstrikes <- rep(0,nrow(standardized_deviations))
# strikes for unequivocally worst: 1
newstrikes <- newstrikes + c(
    standardized_deviations$par1>0.9*max(standardized_deviations$par1) & 
    standardized_deviations$par2>0.9*max(standardized_deviations$par2) & 
    standardized_deviations$par3>0.9*max(standardized_deviations$par3) & 
    standardized_deviations$par4>0.9*max(standardized_deviations$par4))
# strikes for in no dimension best: 1
newstrikes <- newstrikes + c(
    standardized_deviations$par1>1.1*min(standardized_deviations$par1) & 
    standardized_deviations$par2>1.1*min(standardized_deviations$par2) &
    standardized_deviations$par3>1.1*min(standardized_deviations$par3) &
    standardized_deviations$par4>1.1*min(standardized_deviations$par4))
# delete one strike for having the smallest overall (rowmean) standardized deviation
newstrikes <- newstrikes - c(
    rowMeans(standardized_deviations) < 1.1 * min(rowMeans(standardized_deviations)))
# add one strike for having the largest standardized deviation
newstrikes <- newstrikes + c(
    rowMeans(standardized_deviations) > 0.9 * max(rowMeans(standardized_deviations)))
# add one for having an overall standardized deviation larger than 0.5
newstrikes <- newstrikes + c(
    rowMeans(standardized_deviations) > 0.5)
# add one for having any standardized deviation larger than 1
newstrikes <- newstrikes + c(
    apply(standardized_deviations,1,max) > 1)
# add one for having any standardized deviation larger than 2
newstrikes <- newstrikes + c(
    apply(standardized_deviations,1,max) > 2)
# add one for having wrong sign for having the wrong sign for par2
newstrikes <- newstrikes + c(sign(standardized_deviations$par2) != sign(true_th[2]))
# delete one for having no stanrdized deviation larger than 0.2
newstrikes <- newstrikes - c(
    apply(standardized_deviations,1,max) < 0.2)
# delete one for having no stanrdized deviation larger than 0.5
newstrikes <- newstrikes - c(
    apply(standardized_deviations,1,max) < 0.5)
# delete one for having no stanrdized deviation larger than 1
newstrikes <- newstrikes - c(
    apply(standardized_deviations,1,max) < 1)
# delete one for having the best par4
newstrikes <- newstrikes - c(
    standardized_deviations$par4 < 0.9*min(standardized_deviations$par4))
# delete one for having a par4 smaller than 0.5
newstrikes <- newstrikes - c(
    standardized_deviations$par4 < 0.5)
# delete one for having a par4 smaller than 1
newstrikes <- newstrikes - c(
    standardized_deviations$par4 < 1)




print(newstrikes)


igjkh



plot_over_x <- function(which,from,to,par, FUN, ...){
  x <- seq(from,to,length.out=50)
  y <- rep(NA,length(x))
  for (i in 1:length(x)){
    par[which] <- x[i]
    y[i] <- FUN(par,...)
    plot(x,y)
  }
}
plot_from_to <- function(from,to,steps=10, FUN, ...,extend=0){
  sequ <- seq(0-extend,1+extend,length.out=steps)
  y <- rep(NA,length(sequ))
  for (i in 1:length(sequ)){
    par <- from*(1-sequ[i]) + to*sequ[i]
    y[i] <- FUN(par,...)
    plot(sequ,y)
  }
}
plot_over_x(3,-6,2,true_th,FUN=target_function,vdata=vdata,vcv,  keep=keep,prec=2)
vfdsvdfv

results<-rbind(results,
               c(run_1000_new3$par,
                 run_1000_new3$val,
                 target_function(theta=true_th,
                                 vdata=vdata,
                                 vcv=diag(length(keep)),  keep=keep,
                                 prec=10),
                 as.numeric((Sys.time() - ptm),unit="mins"),rseed))

ptm <- Sys.time()

#3 was the same as this:


colnames(results)<-c("p1","p2","p3","p4","fit","time","round")

print(results )
saveRDS(results,"a")

