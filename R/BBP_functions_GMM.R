library(Rfast)
library(stats) #for MLE
library(igraph)
library(foreach)
library(tictoc)
  library(doParallel)
  logistic<-function(x) {
    1/(1+exp(-x))
  }
simulate_BBP<-function(n,delta0,delta1,sigma,distance,kinship,capacity,income,errors=NULL,noiseseed=1,reps=2,parallel=FALSE,computeR=FALSE,plotthis=FALSE) {
  oldseed <- .Random.seed
  #cat("repet:",reps,"\n")
  if (!parallel) {
    finalMatrix<-NULL
    
    for (i in 1:reps) {
      #cat(".",delta0,delta1,sigma,"\n")
      set.seed(noiseseed+i)
      error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
      error <- upper_tri.assign(error,rnorm(n*(n-1)/2,sd=sigma))#make symmetric
      error <- lower_tri.assign(error,lower_tri(t(error)))
      altruism <- 1/(1+exp(-(delta0+delta1*kinship+error)))
      diag(altruism)<-1
      if(mean(upper_tri(altruism)>0.999)>0.5) cat("b")
      eq<-equilibrate_and_plot(altruism=altruism,income=income,modmode=21,capacity=capacity,computeR=computeR,computeCPP=!computeR,plotthis = plotthis)
      ret<-compute_moments(1*(eq$transfers>0),kinship,distance,income,theta=c(delta0,delta1,sigma))
      finalMatrix<-rbind(finalMatrix,ret)
    }
  }
  
  .Random.seed<-oldseed
  rnorm(1) #this should update the seed?
  return(finalMatrix)
}


drawfortheta<-function(theta,kinship,income,distance,...,modelplot=FALSE) {
  n<-length(income)
  error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
  error <- upper_tri.assign(error,rnorm(n*(n-1)/2,sd=exp(theta[3])))#make symmetric
  error <- lower_tri.assign(error,lower_tri(t(error)))
  altruism <- 1/(1+exp(-(theta[1]+theta[2]*kinship+error)))
  diag(altruism)<-1
  if(any(is.nan(altruism))) browser()
  #########emprical vcv####
  eq<-equilibrate_and_plot(altruism=altruism,capacity=-log(theta[4]),income=income,modmode=21,plotthis = modelplot)
  gg<-graph_from_adjacency_matrix((eq$transfers>0)*1)
  plot(gg,...,main="simulated")
  return((eq$transfers>0)*1)
}
gini_from_theta<-function(theta,kinship,income,distance,...,flataltruism=FALSE) {
  n<-length(income)
  error <- matrix(0,nrow=n,ncol=n) #for now, "altruism" is Normal, which is not ideal, given that it is supposed to be in [0,1]
  error <- upper_tri.assign(error,rnorm(n*(n-1)/2,sd=exp(theta[3])))#make symmetric
  error <- lower_tri.assign(error,lower_tri(t(error)))
  altruism <- 1/(1+exp(-(theta[1]+theta[2]*kinship+error)))
  if (flataltruism) {
    altruism[!diag(nrow(kinship))] <- mean(altruism[!diag(nrow(kinship))])
  }
  diag(altruism)<-1
  if(any(is.nan(altruism))) browser()
  #########emprical vcv####
  eq<-equilibrate_and_plot(altruism=altruism,capacity=-log(theta[4]),income=income,modmode=21,plotthis = FALSE)
  return(Rfast::ginis(as.matrix(income+colSums(eq$transfers)-rowSums(eq$transfers))))
}


compute_moments<-function(btransfers,kinship,distance,income,theta) {
  offdiag<-!(diag(nrow(btransfers)))
  g<-igraph::graph_from_adjacency_matrix(btransfers)
  undir<-btransfers+t(btransfers)
  pathlenghts<-igraph::average.path.length(g,directed=FALSE,unconnected=FALSE)/nrow(btransfers)
  fb2<-forestness(undir)
  ib<-intermediation(btransfers)
  sa<-support_fast2(btransfers)
  ra<-recip(btransfers)
  c1<-cor(c(undir[offdiag]),c(kinship[offdiag]))
  c2<-cor(c(undir[offdiag]),c(distance[offdiag]))
  density_<-mean(btransfers[offdiag])
  inc <- matrix(income,nrow=nrow(btransfers),ncol=nrow(btransfers))
  
  
  dat=data.frame(a=c(log(inc/t(inc))[offdiag]),b=sign(log(inc/t(inc)))[offdiag],c=c(kinship[offdiag]))  
  dat$aa<-dat$a*dat$a
  dat$ab<-dat$a*dat$b
  dat$ac<-dat$a*dat$c
  dat$bc<-dat$b*dat$c
  r3<-summary(lm(btransfers[offdiag] ~ ., dat))$sigma^2
  
  da2t=data.frame(a=abs(c(log(inc/t(inc))[offdiag])),b=c(kinship[offdiag]))  

  r4<-summary(lm(undir[offdiag] ~ ., da2t))$sigma^2
  
  c3<-mean(abs(log(inc/t(inc)))[undir==1])
  #sharetransferstoricher<- mean((inc<t(inc))[btransfers>0])
  equated_rest<- btransfers-log(inc/t(inc))-theta[1]-theta[2]*kinship #the extend to which income and kin explain transfers
  sqresidual_proxy<-mean(equated_rest[inc>t(inc)]^2) #only for those where income can explain trnsfers
  #equated_rest2<- btransfers-log(inc/t(inc)+1)-theta[1]-theta[2]*kinship #the extend to which income and kin explain transfers
  #sqresidual_proxy2<-mean(equated_rest2[offdiag]^2) #only for those where income can explain trnsfers
  sqresidual_proxy3<-mean(equated_rest[offdiag]^2) #only for those where income can explain trnsfers
  warning("THIS SHOULD ONLY BE USED IN DEBUG MODE")
  return(cbind(density_,#density                                            X X     alpha                                            
               fb2,#forestness                                              X X     kappa                                
               ib,#intermeidation                                             X     ?????
               sa,#support                                                  X X     kappa                                        
               ra,#reciprocity                                                                                                
               pathlenghts,#pathlenghts                                     X       ????    
               c1,#correatlion between kin and transfers                    X X     beta                      
               c2,#correlation between distance and transfers                                                
               r3,#residual o incometranfers regressed on all                 X     sigma                                 
               r4,#                                           
               c3,           #                                              X       alpha                               
               sqresidual_proxy3#                                           X       sigma    
               #               sharetransferstoricher
               )) #1,1,0,1,0,1,1,0,0,0,1,1
  
}
gg<<-0
g<-function(th,transfers,kinship,distance,...,prec,maxrounds=500){
  #cat("g(",th,")=")
  ret<-c(moment_distance(th=th,transfers=transfers,kinship=kinship,distance=distance,...,maxrounds=maxrounds,prec=prec))
  if (length(ret)==0) { return(Inf)}
  #cat(ret,"\n")
  return(ret)
}

g_new<-function(th,transfers,kinship,distance,...,prec,maxrounds=500){
  #cat("gn(",th,")=")
  ret<-c(moment_distance_new(th=th,transfers=transfers,kinship=kinship,distance=distance,...,maxrounds=maxrounds,prec=prec))
  if (length(ret)==0) { return(Inf)}
  #cat(ret,"\n")
  return(ret)
}
g_dist<-function(th,transfers,kinship,distance,income,vcv,prec=1000,maxrounds=500,noiseseed=1){
  x<-compute_moments_cpp(1*(transfers>0),kinship,distance,income,theta=th)
  
  simx<-simulate_BBP_cpp_parallel(nrow(kinship),th[1],th[2],exp(th[3]),
                                  distance,kinship,matrix(-log(th[4]),nrow(kinship),nrow(kinship)),income,th,reps=prec,noiseseed,maxrounds)
  diff<-sweep(simx,2,colmeans(simx))
  x<-x-colmeans(simx)
 
  x<-x[keep]
  diff<-diff[,keep]
  vcv<-vcv[keep,keep]

  W<-solve(vcv)
  ret<-NULL
  for (i in 1:nrow(diff)) {
    ret<-c(ret,diff[i,]%*%W%*%diff[i,])
  }
  rx<-x%*%W%*%x

  return(mean(c(rx)<ret))
}
moment_distance <- function(th,transfers,kinship,distance,income,prec,noiseseed=1,maxrounds=500,verbose=FALSE,vcv=NULL,keep) {
  inc <- matrix(income,nrow=nrow(kinship),ncol=nrow(kinship))
  offdiag<-!(diag(nrow(kinship)))
  
  
  #x<-compute_moments(1*(transfers>0),kinship,distance,income,theta=th)
  x<-tryCatch(compute_moments_cpp(1*(transfers>0),kinship,distance,income,theta=th), error=function(cond) {return(NA)})
  if (any(is.na(x))) browser()
  
  simx<-simulate_BBP_cpp_parallel(nrow(kinship),th[1],th[2],exp(th[3]),
                                  distance,kinship,matrix(-log(th[4]),nrow(kinship),nrow(kinship)),income,th,prec,noiseseed,maxrounds)
  diff<-sweep(simx,2,x)
  if (is.null(vcv)) {
    vcv<-var(diff)
    diag(vcv)[which(diag(vcv)==0)]<-0.00000001 #replace 0 diagonal elements
    cat("approximating VCV via simulations\n")
  }
  if (verbose) print(rbind(t(x),colmeans(simx),keep,colmeans(diff)))  
  diff<-diff[,keep]
  vcv<-vcv[keep,keep]
  
  ret<-tryCatch({Rfast::colmeans(diff)%*%solve(vcv)%*%Rfast::colmeans(diff)},error=function(cond) {return(Inf)})
  if (is.null(ret)) browser()
  if (is.na(ret)) browser()
  if (ret==Inf) browser()
  return(ret)
}

moment_distance_new <- function(th,transfers,kinship,distance,income,prec,noiseseed=1,maxrounds=500,verbose=FALSE,vcv=NULL,
                                keep
) {
  x<-compute_moments(1*(transfers>0),kinship,distance,income,theta=th)
  
  
  simx<-simulate_BBP_cpp_parallel(nrow(kinship),th[1],th[2],exp(th[3]),
                                  distance,kinship,matrix(-log(th[4]),nrow(kinship),nrow(kinship)),income,th,reps=prec,noiseseed,maxrounds)
  diff<-sweep(simx,2,x)
  if (is.null(vcv)) {
    vcv<-var(diff)
    diag(vcv)[which(diag(vcv)==0)]<-0.00000001 #replace 0 diagonal elements
    cat("approximating VCV via simulations\n")
  }  
  diff<-diff[,keep]
  vcv<-vcv[keep,keep]
  W<-solve(vcv)
  ret<-NULL
  for (i in 1:nrow(diff)) {
    ret<-c(ret,diff[i,]%*%W%*%diff[i,])
  }
  return(mean(ret))
}
lognormal<-function(n,mean,sd) {
  exp(rnorm(n)*sd+mean) - exp(mean+sd^2/2)
}


#utlity function (vector valued)
utlity <- function(consumption,altruism,verbose=FALSE) {
  options(warn=-1)
  cutility=log(consumption)
  cutility[is.nan(cutility)]<- -9999999 #replace -Inf/NaN, because it causes annoyingwarnings
  options(warn=0)
  if(verbose) {print(cutility)}
  utlity=cutility%*%altruism
  
  return(utlity)
}

equilibrate_analytically <- function(altruism,income,capacity,starttransfers=NULL) {
  se = FALSE
  if (is.null(starttransfers)) {
    transfers <- matrix(0,nrow=nrow(altruism),ncol=nrow(altruism))
  } else {
    transfers<-starttransfers
  }
  
  for (r in 0:2000) {
    
    updates<-0
    previoustransfers<-transfers
    #cat("R:")
    for (i in 1:n) { #for each player
      #find the best response
      BR<-BBP_get_BR_analytically_smarter(i=i,transfers=transfers,income=income,altruism = altruism, capacities = capacity)
      
      #cat(i, "(",BR,")")
      #update transfers, only keep net transfers
      if (max(abs(transfers[i,] - BR)) >= 0.00001 | any((transfers[i,]>0) != (BR>0))) {
        #cat("updated,")00001
        updates=updates+1;
        if (se) {browser()}
      }
      transfers[i,] <- BR
      #net out transfers
      transfers<-pmax(transfers-t(transfers),0)
      
    }
    
    if (updates==0) {break}
  }
  if (updates>0) {cat("Best responses did not converge to a NE, probably you need to increase the rounds."); return(FALSE)} #else {cat("Stopped after ",r," rounds. Found a/the nash equilibium\n")}
  return(transfers)
}
#this is required for the optimizer, finding the best response. returns the (negative) utility for individual i
equilibrate <- function(altruism,income,capacity,starttransfers=NULL) {
  se = FALSE
  if (is.null(starttransfers)) {
    transfers <- matrix(0,nrow=nrow(altruism),ncol=nrow(altruism))
  } else {
    transfers<-starttransfers
  }
  
  for (r in 0:2000) {
    
    updates<-0
    previoustransfers<-transfers
    #cat("R:")
    for (i in 1:n) { #for each player
      #find the best response
      BR<-BBP_get_BR(i=i,transfers=transfers,income=income,altruism = altruism, capacities = capacity)
      #cat(i, "(",BR,")")
      #update transfers, only keep net transfers
      if (max(abs(transfers[i,] - BR)) >= 0.00001 | any((transfers[i,]>0) != (BR>0))) {
        #cat("updated,")00001
        updates=updates+1;
        if (se) {browser()}
      }
      transfers[i,] <- BR
      #net out transfers
      transfers<-pmax(transfers-t(transfers),0)
      #test if BRs are sufficiently stable over the round
    }
    
    if (r%%20==0|updates==0|TRUE) {cat("Round",r,"... (",updates," nodes updated their transactions)\n")}
    if (updates==0) {break}
  }
  if (updates>0) {cat("Best responses did not converge to a NE, probably you need to increase the rounds."); return(FALSE)} #else {cat("Stopped after ",r," rounds. Found a/the nash equilibium\n")}
  return(transfers)
}
equilibrate_and_plot<-function(altruism,income,seed=NULL,subtitle=NULL,coords=NULL,capacity=Inf,plotthis=FALSE,modmode=21,computeR=FALSE,computeCPP=TRUE) {
  n<-nrow(altruism)
  if(mean(upper_tri(altruism)>0.999)>0.5) cat("a")
  if (computeCPP) {
    transfers<-equilibrate_cpp_fast8_smarter(altruism,income,matrix(1,nrow(altruism),ncol(altruism))*capacity)
  }
  
  if (computeR) {
    transfers_R<-equilibrate(altruism,income,matrix(1,nrow(altruism),ncol(altruism))*capacity)
    #browser()
    if (!computeCPP){
      transfers <- transfers_R
    }
  }
  if (computeR&computeCPP) {
    if (any(abs(transfers-transfers_R)>0.001) | any((transfers>0.0001 )!=(transfers_R>0.001))) browser() else cat("all good\n");
  }
  
  if (length(c(transfers))==1) browser()
  #if (any(abs(transfers-transfers_cpp)>0.0001) | any((transfers>0)!=(transfers_cpp>0)))   browser()
  ##################### Visualization #####################       2.3963   0.5893   6.4229   8.7627   7.7891   7.9731   4.5527   4.1008   8.1087   6.0493 
  if (plotthis) {
    consumption = income+colSums(transfers)-rowSums(transfers)
    #Plot the altruism network
    g<-graph_from_adjacency_matrix(altruism-diag(n),weighted=TRUE)
    E(g)$width <- E(g)$weight*4 + 1 # offset=1
    c_scale <- colorRamp(c('white','green'))
    E(g)$color = apply(c_scale(E(g)$weight^2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    vertex_attr(g, "label") <- round(consumption,1)
    #E(g)$label<-round(E(g)$weight,2)
    
    #prepare to plot the transactions
    gt<-graph_from_adjacency_matrix(transfers*c(1,-1)[1+(transfers==capacity)],weighted=TRUE) #negative means its at the capacity limit
    vertex_attr(gt, "label") <- paste0(V(gt),": ",round(consumption,1))
    if (is.null(coords)) {
      coords <- layout_nicely(graph_from_adjacency_matrix(altruism-diag(n)+5*(transfers>0),weighted=TRUE)) # layout_with_fr layout_with_drl
    }
    E(gt)$color<-c(rgb(0,0,0),rgb(1,0,0))[1+(E(gt)$weight<0)]
    #output
    V(g)$label.cex<-V(gt)$label.cex <- 0.8
    colfunc <- colorRampPalette(c("white", "red"))
    V(gt)$color  <-V(g)$color  <- colfunc(10)[floor((income-min(income))/(max(income)-min(income))*9)+1]
    plot(g, layout = coords, edge.arrow.size=0, vertex.size = 11,sub=subtitle)
    legendpoints = quantile(income,seq(0, 1, .5))
    legend('topleft',legend=round(legendpoints,1),fill=colfunc(10)[floor((legendpoints-min(income))/(max(income)-min(income))*9)+1],title="Cons. indicated inside nodes\n\nPre-transfer income")
    Sys.sleep(0)
    E(gt)$label<-round(E(gt)$weight,1)
    E(gt)$label.cex<-0.6
    plot(gt, layout = coords,add=TRUE, vertex.size = 11 )
  }
  return(list(transfers=transfers,coords=coords))
}
consumption_weights<-function(alphas,transferdirections,t_conmat) {
  
  
  a<-consumption_weights_cpp(alphas,transferdirections,t_conmat)
  #b<-consumption_weights_old(alphas,transferdirections,t_conmat)$c
  #if (max(abs(a-b))>1e-10) browser()
  return(c(a))
}

consumption_weights_old<-function(alphas,transferdirections,t_conmat) {
  #    alphas=altruism
  n<-nrow(transferdirections)
  transferdirections<-transferdirections-t(transferdirections)
  transferdirections[transferdirections==0]<-NA
  alphas<-pmin(pmax(alphas,1e-10),1)
  care<-(alphas)^transferdirections
  care[is.na(care)]<-0
  care0<-care
  
  
  cw<-matrix(rep(NA,n*n),nrow=n,ncol=n)
  
  for(zz in 1:n) {
    l<-is.na(cw)
    if (all(0==(care[l]))) {break} #there's nothing to get anymore
    ####browser()
    cw[l]<-care[l]
    cw[cw==0]<-NA
    care<-care%*%care0
    diag(care)<-0
    
    
  }
  
  diag(cw)<-1 #CW tells us about how consumptions must relate to each other to fit the equality constraint for linked hhs
  cw[!t_conmat]<-NA
  consumption_fractions<-c(colMeans(cw/rowMeans(cw,na.rm=T),na.rm=T))
  
  #normalize within components, so that average cw = 1 and return as vector
  return(list(c=consumption_fractions,runs=zz)) #colmeans might be non-performant and could be replaced 
}
BBP_T_from_atY_plain<-function(alphas,transferstructure,incomes) {
  #cat("a")
  #a<-BBP_T_from_atY_plain_old(alphas,transferstructure,incomes)
  # cat("c")
  b<-BBP_T_from_atY_plain_cpp(alphas,1*(transferstructure>0),incomes)
  #  cat("b")
  #if (any(a!=b)) browser()
  return(b)
}
BBP_T_from_atY_plain_old<-function(alphas,transferstructure,incomes) {
  transferstructure<-1*(transferstructure>0)
  t_components<-components(graph_from_adjacency_matrix(transferstructure>0))
  t_components_matrix<-matrix(0,nrow=n,ncol=t_components$no)
  t_components_matrix[cbind(1:n,t_components$membership)]<-1
  t_components_csize<-t_components$csize
  t_conmat <- matrix(0,nrow(alphas),nrow(alphas))
  for(cc in 1:t_components$no) {
    t_conmat[t_components$membership==cc,t_components$membership==cc]<-1
  }
  
  consumptions <- BBP_c_from_atY_cpp(alphas,transferstructure,c(incomes),t_components_matrix,t_components_csize,t_conmat)
  #print(consumptions)
  Tr<-BBP_T_from_tYc_cpp(transferstructure,incomes,consumptions,matrix(0,length(incomes),length(incomes)))
  return(Tr)
}
BBP_TC_from_atY<-function(alphas,transferstructure,incomes,equilibrium_check=FALSE,t_components=NULL){
  transferstructure<-1*(transferstructure>0)
  if (is.null(t_components)) {
    t_components<-components(graph_from_adjacency_matrix(transferstructure>0))
  }
  t_components_matrix<-matrix(0,nrow=n,ncol=t_components$no)
  t_components_matrix[cbind(1:n,t_components$membership)]<-1
  t_components_csize<-t_components$csize
  t_conmat <- matrix(0,nrow(alphas),nrow(alphas))
  for(cc in 1:t_components$no) {
    t_conmat[t_components$membership==cc,t_components$membership==cc]<-1
  }
  
  consumptions <- BBP_c_from_atY(alphas,transferstructure,incomes,t_components_matrix,t_components_csize,t_conmat)
  Tr<-BBP_T_from_tYc(transferstructure,incomes,consumptions)
  equilibrium<-NA
  if(equilibrium_check) {
    equilibrium<-BBP_in_equilibrium_YaT(transfer=Tr,income=incomes,altruism=alphas) 
  }
  return(list(
    consumption=consumptions,
    transfers=Tr,
    equilibrium=equilibrium
  ))
}
Rcomponents<-function(x) {
  return(igraph::components(graph_from_adjacency_matrix(x)))
}
BBP_c_from_atY <- function(alphas,transferstructure,incomes,t_components_matrix, t_components_csize,t_conmat) {
  a<-BBP_c_from_atY_old(alphas,transferstructure,incomes,t_components_matrix, t_components_csize,t_conmat)
  b<-BBP_c_from_atY_cpp(alphas,transferstructure,incomes,t_components_matrix, t_components_csize,t_conmat)
  if (any(a!=b)) browser()
  return(b)
}
BBP_c_from_atY_old <- function(alphas,transferstructure,incomes,t_components_matrix, t_components_csize,t_conmat) {
  #get component-wise average incomes
  component_incomes<-c(incomes)%*%t_components_matrix / t_components_csize
  
  #get within-component distribution from alphas (only works via tree structure)
  cw<-consumption_weights(alphas,transferdirections=transferstructure,t_conmat=t_conmat) #components could be returned as a by-procut here
  
  return(c(component_incomes%*%t(t_components_matrix)*cw))
}


BBP_in_equilibrium_YaT <-function(transfers,income,altruism,capacities=Inf,tolerance=1e-4) {
  #a<-BBP_in_equilibrium_YaT_R(transfers,income,altruism,capacities,tolerance)
  #cat("zz")
  b<-BBP_in_equilibrium_YaT_cpp(transfers,income,altruism,matrix(1,nrow(transfers),ncol(transfers))*capacities,tolerance)
  #if (a!=b) browser();
  return(b)
}
BBP_in_equilibrium_YaT_R <-function(transfers,income,altruism,capacities=Inf,tolerance=1e-4) {
  
  consumption = income+colSums(transfers)-rowSums(transfers)
  con <- matrix(consumption,nrow=nrow(transfers),ncol=nrow(transfers))
  offdiag<-!diag(nrow(transfers))
  
  eqcondition <- (abs(t(con)/(con)-altruism)<tolerance & transfers>0) | (transfers==0 & t(con)/(con)>=altruism)
  if (any(is.na(eqcondition))) browser()
  if (all(eqcondition[offdiag])) {
    return(TRUE)
  }else{
    return(FALSE)
  }  
}

mynegutility <- function(mytransfers,i,transfers,altruism,income) {
  
  return(mynegutility_cpp(t(mytransfers),i,transfers,altruism,income))
}
mynegutility_old <- function(mytransfers,i,transfers,altruism,income) {
  if (any(mytransfers<0))
    return(9999)
  transfers[i,] <-  mytransfers
  consumption <- income+colSums(transfers)-rowSums(transfers)
  return(-utlity(consumption,altruism)[i])
}
#foobar <-0
BBP_get_BR <- function(i,transfers,income,altruism,capacities=99999) {
  consumption = income+colSums(transfers)-rowSums(transfers)
  #BR2<-optim(par=transfers[i,],fn=mynegutility_old,i=i,transfers=transfers,income=income,altruism=altruism,method="L-BFGS-B",lower=rep(0,n), upper=pmin(consumption[i],capacities))$par
  #BR3<-optim(par=transfers[i,],fn=mynegutility_old,gr=mynegutility_gr,i=i,transfers=transfers,income=income,altruism=altruism,method="L-BFGS-B",lower=rep(0,n), upper=pmin(consumption[i],capacities))$par
  BRR<-nlminb(start=transfers[i,],objective=mynegutility_old,i=i,transfers=transfers,income=income,altruism=altruism,lower=rep(0,n), upper=pmin(consumption[i],capacities))$par
  #print(BRR)
  #BRRa<-BBP_get_BR_analytically(i,transfers,income,altruism,matrix(1,nrow(transfers),ncol(transfers))*capacities)
  #BRRs<-BBP_get_BR_analytically_smarter(i,transfers,income,altruism,matrix(1,nrow(transfers),ncol(transfers))*capacities)
  BRR[i]<-0
  #if (max(abs(BRR-BRRs))>4e-3|any((BRR>0)!=(BRRs>0))) {print(rbind(BRR,BRRs));browser()}
  #BRC<-BBP_get_BR_analytically_cpp(i,transfers,income,altruism,matrix(1,nrow(transfers),ncol(transfers))*capacities)
  #foobar <<- foobar-1
  #if (foobar<0) browser()
  #transfers[i,]<-BRC
  #consumption = income+colSums(transfers)-rowSums(transfers)
  #cat(i,"updates giving to 29", transfers[i,29], " now 29 has",consumption[29],"\n")
  if(is.nan(mynegutility(BRR,i,transfers,altruism,income))) browser()
  #only do it if it's strictly better:
  if (mynegutility(BRR,i,transfers,altruism,income)>mynegutility(transfers[i,],i,transfers,altruism,income)) {
    BRR<-transfers[i,]
  }
  return(c(BRR))
}
BBP_get_BR_analytically <- function(i,transfers,income,altruism,capacities=99999) {
  mytransfers<-rep(0,ncol(transfers))
  transfers[i,]<-mytransfers
  consumption <- income+colSums(transfers)-rowSums(transfers)
  marginal_util<-(altruism[i,]/pmax(consumption,0.00000001))
  corder<-setdiff(order(marginal_util,decreasing =TRUE),i)
  cat("\n\nI will at most give to",which(marginal_util>1/consumption[i]),"\n")
  #cat("I,",i,", want to give to people in the folowing order", consumption ,"\n")
  for(giveto in 1:(nrow(transfers)-1)) {
    cat("To include ",corder[giveto],"marginal util after giving have to be at most",marginal_util[corder[giveto]], "\n")
    whotogiveto = corder[1:giveto] 
    cat("i still need to optimize for ", whotogiveto,"\n") #can we find those before computing everything?
    CCC=sum(consumption[whotogiveto])
    AAA=sum(altruism[i,whotogiveto])
    TTT=(CCC-AAA*consumption[i])/(-1-AAA)
    
    suggestedtransfers<-altruism[i,whotogiveto]*(consumption[i]-TTT)-consumption[whotogiveto]
    cat("if i could I would give")
    print(rbind(whotogiveto,suggestedtransfers))
    if (any(suggestedtransfers<0)) { #that means we found the full set of peope to give to, but what's the optimum if we want to maintain non-negativity?
      break;
    }
    mytransfers[whotogiveto] <- suggestedtransfers
  }
  return(mytransfers)
}

BBP_get_BR_analytically_smarter <- function(i,transfers,income,altruism,capacities=99999) {
  mytransfers<-rep(0,ncol(transfers))
  transfers[i,]<-mytransfers
  consumption <- income+colSums(transfers)-rowSums(transfers)
  marginal_util<-(altruism[i,]/pmax(consumption,0.00000001))
  corder<-setdiff(order(marginal_util,decreasing =TRUE),i)
  #i should give to anyone who has a larger marginal utiltiy than I
  #cat("\n\nI will at most give to",which(marginal_util>1/consumption[i]),"\n")
  #cat("I,",i,", want to give to people in the folowing order", consumption ,"\n")
  for(giveto in 1:(nrow(transfers)-1)) {
    #cat("To include ",corder[giveto],"marginal util after giving have to be at most",marginal_util[corder[giveto]], "\n")
    whotogiveto = corder[1:giveto] 
    maxout<-NULL
    someoneisnotatcorner<-TRUE
    notoptized<-TRUE
    #round<-0
    consumptionaftermaxtransfers<-consumption
    while (someoneisnotatcorner&notoptized)
    {
      #round=round+1
      #if(round>1) print("toac") 
      #cat("i still need to optimize for ", whotogiveto,"\n") #can we find those before computing everything?
      CCC=sum(consumptionaftermaxtransfers[whotogiveto])
      AAA=sum(altruism[i,whotogiveto])
      TTT=(CCC-AAA*consumptionaftermaxtransfers[i])/(-1-AAA)
      suggestedtransfers<-altruism[i,whotogiveto]*(consumptionaftermaxtransfers[i]-TTT)-consumptionaftermaxtransfers[whotogiveto]
      #cat("if i could I would give")
      #print(rbind(whotogiveto,suggestedtransfers,capacities[i,whotogiveto]))
      excessive<-which(suggestedtransfers>capacities[i,whotogiveto])
      newmaxout<-whotogiveto[excessive]
      if (length(newmaxout)==0) {
        notoptized<-FALSE
      } else {
        maxout<-c(maxout,newmaxout)
        whotogiveto<-setdiff(whotogiveto,newmaxout)
        consumptionaftermaxtransfers[i]=consumptionaftermaxtransfers[i]-sum(capacities[i,newmaxout])
        suggestedtransfers<-suggestedtransfers[!excessive] #in case i move on now, those are the new transfers? It could be that this is alwys empty, no?
      }
      cat("maxout=",maxout,"\nwhotogiveto",whotogiveto,"\nsuggestedtransfers",suggestedtransfers,"\nTTT=",TTT,"\nAAA=",AAA,"\nCCC=",CCC,"\n\n")
      if (length(whotogiveto)==0) someoneisnotatcorner<-FALSE
      #cat("this exceeds limits for ", maxout,"\n") #can we find those before computing everything?
    }
    if (any(suggestedtransfers<0)) { #that means we found the full set of peope to give to, but what's the optimum if we want to maintain non-negativity?
      #if(giveto>1) 
      #  cat("I gave to ",whotogiveto[1:(giveto-1)],"a total of",sum(mytransfers),TTT,"marginal util is",1/(consumption[i]-sum(mytransfers)),"\n\n")
      #else
      #  cat("dontgive\n\n")
      break
    }
    mytransfers[c(whotogiveto,maxout)] <- c(suggestedtransfers,capacities[i,maxout])
  }
  return(mytransfers)
}

BBP_T_from_tYc<- function(transferstructure,incomes,consumptions,Tr=NULL) {
  if (is.null(Tr)) {
    Tr<-matrix(0,length(incomes),length(incomes))
  }
  #a<-BBP_T_from_tYc_old(transferstructure,incomes,consumptions,Tr)$transfers
  b<-BBP_T_from_tYc_cpp(transferstructure,incomes,consumptions,Tr)
  #if (any(b!=a)) browser()
  return(b)
}
BBP_T_from_tYc_old <- function(transferstructure,incomes,consumptions,Tr=NULL,depth=0) {
  #probably this can be significantly sped up by refusing to wokr on cases where transfers are inconsistent with consumption
  #find nodes that have only one t
  if (is.null(Tr)) {
    Tr <- matrix(0,nrow=nrow(transferstructure),ncol=nrow(transferstructure))
  }
  #receivers
  #print("a")
  indegree<-colSums(transferstructure)
  outdegree<-rowSums(transferstructure)
  intermediateincomes<-incomes+colSums(Tr)-rowSums(Tr)
  
  canTupdate<-which(indegree!=1 | outdegree!=0)
  netflows_to_receivers<- consumptions-intermediateincomes #otherwise this counts double
  
  
  if(any(indegree>0 & outdegree==0 & netflows_to_receivers<0)) {
    #cat("Rnot an eq\n")
    return(list(transfers=Tr,equilibrium=FALSE))
  }
  
  netflows_to_receivers[canTupdate]<-0
  newTr<-t(t(transferstructure)*netflows_to_receivers)
  
  Tr[newTr>0]<-newTr[newTr>0]
  transferstructure[Tr>0]<-0
  
  #givers
  indegree<-colSums(transferstructure)
  outdegree<-rowSums(transferstructure)
  intermediateincomes<-incomes+colSums(Tr)-rowSums(Tr)
  canTupdate<-which(outdegree!=1 | indegree!=0)
  netflows_from_givers<-consumptions-intermediateincomes #otherwise this counts double
  
  if(any(outdegree>0 & indegree==0 & netflows_from_givers>0)) {
    #cat("Rnot an eqww\n")
    return(list(transfers=Tr,equilibrium=FALSE))
  }
  
  netflows_from_givers[canTupdate]<-0
  newTr<-transferstructure*(-netflows_from_givers)
  Tr[newTr>0]<-newTr[newTr>0]
  transferstructure[Tr>0]<-0
  
  if(max(transferstructure)==0) {
    return(list(transfers=Tr,equilibrium=NA))
  } else if(depth>nrow(transferstructure)) {
    #cat("canthelp to solve this\n")
    return(list(transfers=Tr,equilibrium=FALSE))
  } else {
    return(BBP_T_from_tYc_old(Tr=Tr,transferstructure = transferstructure,incomes,consumptions,depth=depth+1))
  }
}


forestness<-function(adj) {
  if (any(adj!=t(adj))) {
    warning("could it be that you're not passing a symmetric matrix to forestness?")
    browser()
  }
  return(2*(nrow(adj)-component_counts(adj)) / sum(adj))
}
intermediation<-function(adj) {
  radj<-adj*t(adj)
  nradj<-adj-radj
  return(sum(rowMeans(nradj)>0 & colMeans(nradj)>0)/sum(rowMeans(nradj)>0 | colMeans(nradj)>0))
}
recip <-function(adj)
{
  radj<-adj*t(adj)
  return(sum(radj)/sum(adj))
}
support<-function(g) {
  gu<-as.undirected(g)
  pairs_of_friends<-length(E(gu))
  E(gu)$supported<-NA
  for (e in E(gu)) {
    #browser()
    options(warn=-1)
    other_path<-shortest_paths(delete_edges(gu,e),ends(gu,e)[1],ends(gu,e)[2])$vpath
    options(warn=0)
    E(gu)$supported[e]<-any(sapply(other_path,function(path){return(length(path))})==3)
  }
  return(sum(E(gu)$supported)/gsize(gu))
}
support_fast<-function(m) {
  undir<-pmax(m,t(m))
  twopaths = undir%*%undir
  return(mean((twopaths>=1)[undir==1]))
}
support_fast2<-function(m) {
  undir<-pmax(m,t(m))
  twopaths = undir%*%undir
  return(mean(((twopaths[upper.tri(m)])>=1)[undir[upper.tri(m)]==1]))
}
multiplicity<-function(g) {
  pathcount<-matrix(NA,nrow=gorder(g),ncol=gorder(g))
  options(warn=-1)
  for (i in V(g)) {
    for (j in V(g)) {
      pathcount[i,j]<-0
      first_path<-shortest_paths(g,i,j)$vpath[[1]]
      if (length(first_path)==0|i==j) next
      pathcount[i,j]<-1
      for (z in 1:(length(first_path)-1)) {
        other_path<-shortest_paths(delete_edges(g,paste(c(first_path[z],first_path[z+1]),collapse="|")),i,j)$vpath[[1]]
        if (length(other_path)>0) {pathcount[i,j]<-2; break}
      }
    }
  }
  options(warn=0)
  diag(pathcount)<-0
  
  return(sum(pathcount>1)/sum(pathcount>0))
}

# ###########
# fit:
#   density     forest      support     pathlengths correlation1
# [1] 0.03721000 0.84610393 0.05619804 0.21139954 0.15015951
# true:
# [1] 0.03637222 0.83635829 0.02018914 0.24995272 0.15071341
# 
# real moments:
# [1] 0.03777778 0.85294118 0.08823529 0.16421456 0.14973230
# 
# 
# 
# fb2         ib         sa ra pathlenghts        c1           c2
# [1,] 0.03777778 0.8529412 0.06666667 0.08823529  0   0.1642146 0.1497323 -0.001912973
# Browse[1]> x[keep]
