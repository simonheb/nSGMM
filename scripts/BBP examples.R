rm(list=ls())
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(mlrMBO)

setwd("C:/Users/shess/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')
source('R/gridsearch.R')

drawnetexample<-function(delta0_DGP,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,llayout=NULL,showunconstrained=FALSE,prefix=""){
  
  
  altruism <- 1/(1+exp(-(delta0_DGP+delta1_DGP*kinship+error*sigma_DGP)))
  if (showunconstrained) kinship[]<-0
  print(summary(c(altruism)))
  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity_DGP,income=income,modmode=21,plotthis = F)
    equ<-equilibrate_and_plot(altruism=altruism,capacity=99999999,income=income,modmode=21,plotthis = F)
    tnu<-simplify(graph_from_adjacency_matrix(equ$transfers,weighted=TRUE))

  an<-as.undirected(simplify(graph_from_adjacency_matrix(altruism,weighted=TRUE,mode = c("max"))))
  E(an)$width <- E(an)$weight^2*10   # offset=1
  colfunc <- colorRampPalette(c("#00000000", "#000000FF"), alpha = TRUE)
  E(an)$color <- colfunc(100)[round(E(an)$weight*100)+1]

  tn<-simplify(graph_from_adjacency_matrix(eq$transfers,weighted=TRUE))
  kn<-as.undirected(simplify(graph_from_adjacency_matrix(kinship)))
  kn_i<-as.undirected(simplify(graph_from_adjacency_matrix(matrix(0,nrow(kinship),nrow(kinship)))))
  if (is.null(llayout))
    llayout = layout_nicely(as.undirected(simplify(graph_from_adjacency_matrix(altruism+2*kinship+1*eq$transfers,weighted=T))))
  
  
  
  V(tnu)$size  <- V(tn)$size  <-V(kn)$size  <-V(an)$size  <- (income-min(income))/(max(income)-min(income))*20+12
  V(tnu)$color  <-V(tn)$color  <-V(kn)$color  <-V(an)$color  <- "white"
  V(kn_i)$label<-round(income,1)
  V(kn_i)$label.dist<-V(kn)$size/10+0.7
  print(a<-paste0(c("D:/Dropbox/Gambia/working folder/Marcel Allocation/Draft/current - Simon/figures/exampleplot",prefix,delta0_DGP*10,delta1_DGP*10,sigma_DGP*10,10*capacity_DGP,".pdf"),collapse=""))
  pdf(a,width=3.3*(3+showunconstrained),height=4)
  par(mfrow=c(1,3+showunconstrained))
  if (showunconstrained) {
    plot(kn,layout=llayout,edge.width=4,main="Income")
  } else {    
    plot(kn,layout=llayout,edge.width=4,main="Kin and Income")
  }
  plot(kn_i,layout=llayout,add=TRUE, label.dist=V(kn_i)$label.dist, vertex.color = NA, vertex.frame.color=NA)

  plot(an,layout=llayout,main="Altruism")
  if (is.numeric(E(tn)$weight))
    E(tn)$label<-round(E(tn)$weight,1)
   E(tnu)$label<-round(E(tnu)$weight,1)
  if (showunconstrained) {
    plot(tnu,layout=llayout,main="Transfers w/o constraints" ,edge.arrow.size=0.9)
    plot(tn,layout=llayout,main="Transfers w/ constraints",edge.arrow.size=0.9)
  } else {
    plot(tn,layout=llayout,main="Transfers",edge.arrow.size=0.9)
  }
  dev.off()
  return(llayout)
}

seed<-1#round(runif(1)*100)
set.seed(seed)
############# Set Simulation Parameters and Draw Random Altruism Network#############
#population size 
n<-10
delta0_DGP<- -1.8
delta1_DGP<- 2

sigma_DGP<-1.3
capacity_DGP<-99
true_th<-c(delta0_DGP,delta1_DGP,log(sigma_DGP),capacity_DGP)

lon<-rnorm(n)
lat<-rnorm(n)
distance<-as.matrix(dist(cbind(lon,lat),"euclidean",TRUE,TRUE)) #distances between coordinates drawn from a bivariate normal
distance<-distance/max(distance)
diag(distance)<-0

kinship <- matrix(rbinom(n*n,1,0.046),nrow=n)
kinship <- lower_tri.assign(kinship,lower_tri(t(kinship))) #make symmetric
diag(kinship)<-1
  
error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
error <- upper_tri.assign(error,rnorm(n*(n-1)/2))#make symmetric
error <- lower_tri.assign(error,lower_tri(t(error)))

orders<-c(6L, 1L, 2L, 3L, 8L, 9L, 5L, 10L, 4L, 7L)
distance<-distance[orders,orders]
error<-error[orders,orders]
kinship<-kinship[orders,orders]
income <-c(4.906, #1
           5.489, #2
           5.907, 
           6.149, #4
           6.556, #5
           6.787, 
           17.809, #7
           6.348, 
           13.321, 
           31.141)
error<-  structure(c(0, -0.57, -0.14, 2.18, -1.28, 0.16, 1.69, -0.46, -1.52, -0.54,
                     -0.57, 0, 0.48, -0.71, 1.21, -1.22, -0.44, 0.72, -0.93, 0.33,
                     -0.14, 0.48, 0, 0.61, 1.16, -0.47, 0, 0.91, -1.25,  1.06,
                     2.18, -0.71, 0.61, 0, 0.7, -0.62, 0.07, 0.38, 0.29, -0.3,
                     -1.28, 1.21, 1.16, 0.7, 0, 1.77, 0.56, -0.65, 1.59, -0.57,
                     0.16, -1.22, -0.47, -0.62, 1.77, 0, -0.91, -0.21, 0.04, -0.65,
                     1.69, -0.44, 0, 0.07, 0.56, -0.91, 0, -0.64, -0.59, 0.27,
                     -0.46, 0.72, 0.91, 0.38, -0.65, -0.21, -0.64, 0, 1.68, 1.43,
                     -1.52, -0.93, -1.25, 0.29, 1.59, 0.04, -0.59, 1.68, 0, 0.37,
                     -0.54, 0.33, 1.06, -0.3, -0.57, -0.65, 0.27, 1.43, 0.37, 0), .Dim = c(10L, 10L))

layout<-structure(c(2.19177352918267, 0.927656391577087, 0.530294198246127, 
                    1.6038303046075, 0.959225024252714, 0.357827028068436, 1.52513782239813, 
                    -0.0169703077756292, 0.367104103919506, -0.02992720103484, 3.32446205693794, 
                    2.66296327404841, 4.26893158720943, 3.85417334884739, 3.21025286934964, 
                    2.54564190302698, 4.62023584293743, 3.13601025313869, 3.62434499631413, 
                    3.86593434362282), .Dim = c(10L, 2L))
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP+2,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP-1,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP+1,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP-1,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP-1,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,1,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,0.5,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,3,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,0.5,kinship,error,income,layout,showunconstrained = TRUE, prefix= "maintext")





  
  
