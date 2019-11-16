rm(list=ls())
library(Rcpp)
library(foreach)
library(doParallel)
library(nSGMM)
library(mlrMBO)

setwd("D:/Dropbox/Gambia/working folder/Marcel Allocation/nSGMM")
source('R/BBP_functions_GMM.R')
source('R/gridsearch.R')

drawnetexample<-function(delta0_DGP,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,llayout=NULL){
  
  
  
  altruism <- 1/(1+exp(-(delta0_DGP+delta1_DGP*kinship+error*sigma_DGP)))
  eq<-equilibrate_and_plot(altruism=altruism,capacity=capacity_DGP,income=income,modmode=21,plotthis = F)
  
  an<-as.undirected(simplify(graph_from_adjacency_matrix(altruism,weighted=TRUE,mode = c("max"))))
  E(an)$width <- E(an)$weight^2*10   # offset=1
  colfunc <- colorRampPalette(c("#00000000", "#000000FF"), alpha = TRUE)
  E(an)$color <- colfunc(100)[round(E(an)$weight*100)+1]

    
  tn<-simplify(graph_from_adjacency_matrix(eq$transfers,weighted=TRUE))
  kn<-as.undirected(simplify(graph_from_adjacency_matrix(kinship)))
  kn_i<-as.undirected(simplify(graph_from_adjacency_matrix(matrix(0,nrow(kinship),nrow(kinship)))))
  if (is.null(llayout))
    llayout = layout_nicely(as.undirected(simplify(graph_from_adjacency_matrix(altruism+2*kinship+1*eq$transfers,weighted=T))))
  
  
  
  V(tn)$size  <-V(kn)$size  <-V(an)$size  <- (income-min(income))/(max(income)-min(income))*20+12
  V(tn)$color  <-V(kn)$color  <-V(an)$color  <- "white"
  V(kn_i)$label<-round(income,1)
  V(kn_i)$label.dist<-V(kn)$size/10+0.7
  #pdf(paste0(c("D:/Dropbox/Gambia/working folder/Marcel Allocation/Draft/figures/exampleplot",delta0_DGP*10,delta1_DGP*10,sigma_DGP*10,10*capacity_DGP,".pdf"),collapse=""), width=10, height=4)
  print(paste0(c("D:/Dropbox/Gambia/working folder/Marcel Allocation/Draft/figures/exampleplot",delta0_DGP*10,delta1_DGP*10,sigma_DGP*10,10*capacity_DGP,".pdf"),collapse=""))
  par(mfrow=c(1,3))
  
  plot(kn,layout=llayout,edge.width=4,main="Kin and Income")
  plot(kn_i,layout=llayout,add=TRUE, label.dist=V(kn_i)$label.dist, vertex.color = NA, vertex.frame.color=NA)
  
  plot(an,layout=llayout,main="Altruism")
  if (is.numeric(E(tn)$weight))
    E(tn)$label<-round(E(tn)$weight,1)
  plot(tn,layout=llayout,main="Transfers")
  #dev.off()
  
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

income = round(exp(rnorm(n)+2),2)

distance<-distance[order(income),order(income)]
error<-error[order(income),order(income)]
kinship<-kinship[order(income),order(income)]
income<-income[order(income)]

layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income)
layout<-drawnetexample(delta0_DGP+2,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP-1,delta1_DGP,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP+1,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP-1,sigma_DGP,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,3,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP-1.2,capacity_DGP,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,1,kinship,error,income,layout)
layout<-drawnetexample(delta0_DGP,delta1_DGP,sigma_DGP,0.5,kinship,error,income,layout)

#altruism<-base$altruism
#transfers<-base$transfers



  
  
