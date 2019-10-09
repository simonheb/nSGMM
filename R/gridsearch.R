colMadss<-function(x) {colMeans(abs(sweep(x,2,colMeans(x))))}
library(scales)
dimensionwiseGridSearch<-function (fun,...,
                                   stepsize=10, #into how many intervals is each parametre split
                                  # stepoverlap=1, #determines how many minima are kept per round
                                   #stepexpand=1,#how big is the interval around each kept minimum that is going to be re-used after
                                   stepdepth=10,#how many rounds are needed?
                                   plotit=FALSE,lower,upper,
                                  stepexpand=1.5,
                                  pruneratio=0.5
                                  ) {
  if (stepsize<4) message("stepsizes smaller than 4 are discouraged")
  if (1-pruneratio<=1/stepsize) message("prunerations that will only leave a single point are discouraged")
  #loop depth times 
  global_lower<-lower
  global_upper<-upper
  for (i in 1:stepdepth){
    cat("round:", i, "\n")
    #loop over parameters
    for (p in 1:length(lower)){
      current<- rowMeans(cbind(lower,upper))
      steps<- (upper[p]-lower[p])*((0:(stepsize-1))/(stepsize-1)*stepexpand-(stepexpand-1)/2 ) +lower[p]
      steps<-unique(pmin(pmax(steps,global_lower[p]),global_upper[p]))
      x<-NULL
      y<-NULL
      for(s in steps) {
        x<-c(x,s)
        x_all<-current
        x_all[p]<-s
        y<-c(y,fun(x_all,prec=i^3,...))
      }
      if (plotit) plot(x,y,main=p)
      keepers<-which(y<=quantile(y,1-pruneratio))
      oldspan<-(upper[p]-lower[p])
      cat("parameter",p,". reduction to", (max(x[keepers])-min(x[keepers]))/oldspan,"\n")
      lower[p]<-ifelse(min(x[keepers])>lower[p],min(x[keepers]), lower[p]-oldspan*stepexpand) #if the max is at the border, expand to that side, otherwise, zoom in
      upper[p]<-ifelse(max(x[keepers])<upper[p],max(x[keepers]), upper[p]+oldspan*stepexpand) #if the max is at the border, expand to that side, otherwise, zoom in
    }
  }
  return(cbind(lower,upper))

}

zoomingGridSearch<-function (fun,...,lower,upper,
                             stepsize=10, #into how many intervals is each parametre split
                             stepoverlap=1, #determines how many minima are kept per round
                             stepexpand=1,#how big is the interval around each kept minimum that is going to be re-used after
                             stepdepth=10,
                             pruneratio=0.5,prunepoints=FALSE,
                             plotit=FALSE) {
  if (any(upper<lower)) browser()
  points = NULL
  nldf<-matrix(c(0,0,0,0),nrow=1)
  p_midpoints = matrix((lower+upper)/2,nrow=1)
  p_span = upper-p_midpoints
  for (step in 1:stepdepth) {
    newlevels<-NULL
    #first rasterize midpoints in order to reduce complexity
    #rounding to avoid calculating similar points twice:
    for (c in 1:ncol(p_midpoints)) {
      p_midpoints[,c]<-plyr::round_any(p_midpoints[,c],(2*p_span[,c]+max(p_midpoints[,c])-min(p_midpoints[,c]))/(2*stepsize))
    } 
    
    for(rows in 1:nrow(p_midpoints)) {
      steps_of_each_parameter<-matrix(((-1*stepsize):(1*stepsize))/(stepsize),ncol=1)%*%matrix(p_span,nrow=1) #I draw 5 times as many, they will be reduced below by the rasterization anyways
      all_possible_combinations_of_deviations<-expand.grid(as.list(as.data.frame(steps_of_each_parameter)))
      all_parameter_combinations<-sweep(all_possible_combinations_of_deviations, 2, p_midpoints[rows,], "+")
      valid_combinations<-all_parameter_combinations[colMeans(t(all_parameter_combinations)>upper | t(all_parameter_combinations)<lower)==0,]
      newlevels<-rbind(newlevels,valid_combinations)
    }

    #this shouldnt be needed, but apparently it is due to numeric imprescision
    for (c in 1:ncol(newlevels)) {
      newlevels[,c]<-plyr::round_any(newlevels[,c],(2*p_span[,c]+max(p_midpoints[,c])-min(p_midpoints[,c]))/(2*stepsize))
    } 
    
    

    if(plotit) {
      plot(rbind(nldf[,1:2],newlevels[,1:2]),col = alpha("white", 0.0), pch=16) 
      points(nldf[,1:2],col = alpha("yellow", 0.4), pch=16) 
      points(newlevels[,1:2],col = alpha("blue", 0.1), pch=16) 
      #points(t(c(9.87654321,1.23456789)),col = alpha("red", 1), pch=16,cex=3) 
      nldf<-rbind(newlevels)
    }

    if(prunepoints) {
      
      newlevels<-plyr::count(newlevels)
      print(table(newlevels$freq))
      cat("cutoff",quantile(newlevels$freq,pruneratio), "\n")
      message(paste0("reducing ", nrow(newlevels)))
      
      newlevels<-newlevels[which(newlevels$freq>=quantile(newlevels$freq,pruneratio)),1:(ncol(newlevels)-1)]

      message(paste0("to ",nrow(newlevels)))
    }
    if(plotit) {
      points(newlevels[,1:2],col = alpha("green", 0.8), pch=19) 
      nldf<-rbind(newlevels)
    }
    newlevels<-unique(as.list(as.data.frame(t(newlevels))))
    for(i in 1:700) {
      if (i%%100==0) cat(i)
      if (i%%100<3) next 
      if (i%%100==99) cat("|")
      cat("-")
      
    }
    #cat(paste(c(rbind(1:7*100, paste(c(replicate(97, "-"),"|"), collapse = ""))),collapse=""))
    p_midpoints = unique(gridSearch(fun,levels=newlevels,prec=step,...,nmin=stepoverlap)$minlevels)
    if(plotit) {
      points(p_midpoints[,1:2],col = alpha("red", 0.8), pch=19) 
      dev.flush()
    }

    p_span<-pmax(p_span/(stepsize)*stepexpand,
                 colMadss(p_midpoints)*stepexpand)
    
    
    cat("found")
    print(p_midpoints)
    #cat("new span\n")
    #print(p_span)
    points<-rbind(points,p_midpoints)
  }
  if(plotit) plot(points[5:step,1],t="l")
  return(p_midpoints[1,])
}
gridSearch<-function (fun, levels, ..., 
          printDetail = TRUE,nmin=1) 
{
  message(paste0("evaluating:",length(levels)))
  results <- lapply(levels, fun, ...)
  results <- unlist(results)
  is <- try(order(results))
  i=is[1:nmin]
  if (inherits(is, "try-error") || any(is.na(is)) || length(i) == 
      0L) {
    warning("cannot compute minimum (NA values in results, ...)")
  }
  else {
    list(minfun = results[i], minlevels = t(matrix(unlist(levels[i]),nrow=length(levels[[1]]))))
  }
}
#f<-function(x,prec=Inf){cat(x,".\n");1.9*(x[1]-0.987654321)*(x[2]-0.123456789)+(x[1]-0.987654321)^2+(x[2]-0.123456789)^2+(x[4]-0.123456789)^2 + (x[3]-0.123456789)^2+rnorm(1)/(prec)}
#gridSearch(f, levels=list(c(1,2),c(2,3),c(9.8765,1.2345)))
#optim(c(0,0,0,0),f, lower=rep(-50,4),upper=rep(50,4),prec=Inf)

#good expeirience with stepexpand=2/3*stepsize
#zoomingGridSearch(f, lower=rep(-50,4),upper=rep(50,4),stepdepth=20,stepsize=3,stepoverlap=8,stepexpand=1.3,plotit=TRUE,prunepoints = TRUE,pruneratio=0.75)
#a<-dimensionwiseGridSearch(f, lower=rep(-50,4),upper=rep(50,4),stepdepth=10,stepsize=12,plotit=TRUE,pruneratio=0.66,stepexpand=1.5)
#print(rowMeans(a))
