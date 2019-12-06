colMadss<-function(x) {colMeans(abs(sweep(x,2,colMeans(x))))}
library(scales)
library(mlrMBO)
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
plot_partial<-function(theta,param=1,minoffset=-1,maxoffset=1,fun,steps=4,...){
  x<- (0:steps)/steps*(maxoffset-minoffset)+minoffset
  print(x)
  xs<-lapply(x,function(x,theta,par){theta[par] = theta[par] +x;theta},theta=theta,par=param)
  plot(x,lapply(xs,fun,...))
}

random_walk_cooling<-function(fun,#the function to be minimized
                              theta,#staring value
                              ...,#other options to be passed to fun
                              stepsize=NULL,#initial stepsize for candidate search 
                              reheatineval=100,#after how many iterations do we consider going back to an old optimum?
                              warmuplength=5,#for how many periods are switching probabilities overproportional to the relative imporvement
                              iter=1000,#how many iterations?
                              momentum=0,#momentum direction (usually not specified)
                              momentumdecay=0.5,#speed at which momentum decays
                              trace=NULL, #initial return value
                              maxiter=0, #where did we start? 
                              tempresetiter=NULL, #when did we reset temperature?
                              minstep=0.01, #mimum stepsize
                              upper,  lower, #limits (important!)
                              precschedule,#a function that returns a integer that is passed as prec to fun
                              true_theta=NULL,#for plotting
                              roundstopartial=200,
                              superverbose=FALSE,
                              current=NULL,currentvalueage=10,
                              lastplottet=Sys.time()-60,
                              bestpast=NULL,
                              true_g=NULL#the value a restet will use
) {
  if(iter==0) {
    
    cat("\n")
    f1<-fun(bestpast$theta,prec=4*precschedule(iter,maxiter),...)
    f2<-fun(theta,prec=4*precschedule(iter,maxiter),...)
    if (f2>f1) {
      cat("returning something I found along the way",bestpast$iter,"\n")
      print(bestpast$theta)
      return(list(theta=bestpast$theta,
                  trace=trace,
                  value=f1,
                  bestpast=bestpast,
                  true_g=true_g))
    }
    print(theta)
    return(list(theta=theta,
                trace=tail(trace,maxiter-bestpast$iter),
                value=f2,
                bestpast=bestpast,
                true_g=true_g))
  }
  
  maxiter=max(iter,maxiter)
  if (is.null(tempresetiter)) tempresetiter=maxiter
  if (is.null(stepsize)) stepsize <-(upper-lower)/100
  if (!any(is.null(true_theta)) & is.null(true_g)) {
    true_g<- fun(th=true_theta,prec=4*precschedule(1,maxiter),...) 
  }  else {
    true_g<- -1
  } 
    
  if(iter%%roundstopartial==50) { #every 100 rounds solve for the partial min
    theta<-partialoptim(fun,theta,prec=precschedule(iter,maxiter),lower=lower,upper=upper,steps=4,maxrounds=precschedule(iter,maxiter)/4,...)
    currentvalueage<-Inf #make sure this gets re-computed.
  }
  
  momentum<-momentum*momentumdecay
  step<-rnorm(length(theta))*stepsize+momentum
  #print(rbind(momentum,stepsize,step))
  newtheta<-theta+step
  theta<-pmax(pmin(theta,upper),lower)
  newtheta<-pmax(pmin(newtheta,upper),lower)
  #step<-newtheta-theta
  
  
  #if
  #old<-max(0.000001,fun(theta,prec=precschedule(iter,maxiter),,...,noiseseed=i))
  if(currentvalueage>=10) { #after recycling the same value 10 times, let's  use a new draw so that we are not getting stuck in a lucky draw from g
    old<-max(0.000001,fun(theta,prec=precschedule(iter,maxiter),noiseseed=iter,...))
    currentvalueage<-0
    cat("=")
  }else {
    old<-current
  }
  #cat(old,"=",current,"\n")
    
  if (old==Inf & maxiter==iter) {browser();stop("starting at a point that evaluates to Inf is unlikely to yield a good result")}
  new<-max(0.000001,fun(newtheta,prec=precschedule(iter,maxiter),...,noiseseed=iter))
  if(old==Inf & new==Inf) {old<-new<-999999999} #because sometimes both are Inf
  if(runif(1)<(min(old/new,1)^((tempresetiter-iter)/warmuplength))) {
    ##sometimes this runs into oblivion again, whenever this happens, get it back home
    if (iter<600 & new>60*max(bestpast$fit,new) | iter<300 & new>20*max(bestpast$fit,new) | new>200*max(bestpast$fit,new)) { #can only be true if bestpas!=NULL
      newtheta<-bestpast$theta
      new<-bestpast$fit
      cat("*")
    }
    effectivestep<-newtheta-theta
    theta<-newtheta
    current<-new
    currentvalueage<-0
    stepsize<-2 * #blow up stepsize to be more probing, 
      pmax(stepsize+
             iter/tempresetiter #add some extra searching early on
           ,minstep)
    if (maxiter-iter>2) stepsize<-pmin(stepsize,(upper-lower)/10) #prevent explosion, e.g. because prolonged warmup
    momentum<-momentum+effectivestep
    cat("|")
    if (new<old | is.null(bestpast))  { #if the new value is better than old
      if (!is.null(bestpast))
        bestpast_fit<-bestpast$fit
      else 
        bestpast_fit<-Inf
      if (new<bestpast_fit) {
        bestpast<-list(theta=newtheta,iter=iter,fit=new)
      } else { #hope is, this would happen rarely  #better than old, but worse than the pastbest, maybe the pastbest was too optimist, so let's re-evaluate it so that we won't be stuck with a lucky draw
        cat("o")
        bestpast$fit<-max(0.000001,mean(c(bestpast$fit,fun(bestpast$theta,prec=precschedule(iter,maxiter),...,noiseseed=iter))))
      }
    }
  } else {
    stepsize<-stepsize/2
    cat(".")
    current<-old
    currentvalueage<-currentvalueage+1
  }
  
  trace=rbind(trace,c(theta,current))
  if(iter%%reheatineval==0) {
    #reset the chain
    if (!is.null(bestpast)) {
      if (bestpast$iter>iter & iter > 50) {
        bestpast$fit<-max(0.000001,fun(bestpast$theta,prec=precschedule(iter,maxiter),...)) #first, update
        if (runif(1)*2<min(1,current/bestpast$fit)) { #if the past is better and with luck, reset.
          cat("[",iter-bestpast$iter,"]") 
          tempresetiter<-iter+warmuplength
          theta<-bestpast$theta
          stepsize<-head(colSd(trace[ceiling(0.9*nrow(trace)):nrow(trace),]),length(theta))
        }
      }
    }
  }
  if((iter-1)%%200==0|iter<=1|as.numeric(Sys.time()-lastplottet,units="mins")>1+(maxiter-iter)/100) {
    if (is.null(true_theta)) {
      plot_trace(trace,c(theta,true_g),bestpast,maxiter)
    } else {
      plot_trace(trace,c(true_theta,true_g),bestpast,maxiter)
    } 
    lastplottet<-Sys.time()
    cat("\n[iter",iter,":c(",paste0(round(theta,3),collapse=","),")=",current,"; pastbest@",maxiter-bestpast$iter,", was",bestpast$fit," ]\n")
    
  }
  return(random_walk_cooling(fun,theta,...,stepsize=stepsize,iter=iter-1,tempresetiter=tempresetiter,momentum=momentum,minstep=minstep,momentumdecay=momentumdecay,trace=trace,maxiter=maxiter,reheatineval=reheatineval,warmuplength=warmuplength,upper=upper,lower=lower,true_theta=true_theta,precschedule=precschedule,bestpast=bestpast,current=current,currentvalueage=currentvalueage,lastplottet=lastplottet,true_g=true_g))
}


partialoptim<-function(fun,#the function to be minimized
                       theta,#staring value
                       ...,#other options to be passed to fun
                       prec=1000,#how many iterations?
                       upper,  lower, #limits (important!)
                       steps=8
) {
  curr<-fun(theta,prec=prec,...)
  if (length(lower)==1) lower = rep(lower,length(theta))
  if (length(upper)==1) upper = rep(upper,length(theta))
  cat("{",round(curr,3))
  loopingpars<-sample(1:length(theta),length(theta))
  backups<-2
  needstobackup<-FALSE
  working<-TRUE
  par_i<-1
  while(working) {
    par<-loopingpars[par_i]
    xs_lower<- c( (1-((1:steps)/steps)^2)*(theta[par]-lower[par])+lower[par])
    xs_upper<- c( ((1:steps)/steps)^2*(upper[par]-theta[par])+theta[par])
    xs<-unique(c(xs_lower,xs_upper))
  #  print(xs)    
    updated<-FALSE
    for (x in xs) {
      tsugg<-theta
      tsugg[par]<-x
      cand<-fun(tsugg,prec=prec,...)
      if(cand<curr){
        theta<-tsugg
        curr<-cand
        updated<-TRUE
      } 
    }
    if (updated) {
      needstobackup<-TRUE
      cat("x")
      par_i<-par_i+1
    } else {
      cat(">")
      par_i<-par_i+1
    }
    if (par_i>length(loopingpars)) {
      if (backups>0 & needstobackup) {
        needstobackup<-FALSE
        backups<-backups-1
        par_i<-1
        cat("|")
      } else {
        working<-FALSE
      }
    }
  }
  cat(round(curr,3),"}")
  return(theta)
}

which.min2 <- function(x, last.index = TRUE, ...){
  if(last.index) max(which(x == min(x, ...,na.rm=TRUE))) else which.min(x)
}

plot_trace<-function(trace,true_theta,bestpast,maxiter) {
  tryCatch({
    par(mfrow=c(ncol(trace),1))
    hideshare<-0.1
    ttrace<-trace
    ttrace[1:ceiling(hideshare*nrow(trace)), 1:(ncol(trace))]<-NA
    #ttrace[1:ceiling((1-(1-hideshare)^3)*nrow(trace)),ncol(trace)]<-NA
    #ttrace[1,1:(ncol(trace)-1)]<-0
    #for (i in 1:(ncol(trace)-1)) {
    #  ttrace[2,i]<-true_theta[i]
    #}
    optloc<-maxiter-bestpast$iter#which.min2(ttrace[,ncol(trace)]);
    ttrace[ttrace[,ncol(trace)]>max(true_theta[5]*2,7*min(ttrace[,ncol(trace)],na.rm=TRUE)),ncol(trace)]<-NA
    plot(ttrace[,1],type="l",ylab="1",     ylim=c( min(0,true_theta[1],tail(ttrace[,1],1),trace[optloc,1])-1,max(0,true_theta[1],tail(ttrace[,1],1),trace[optloc,1])+1)); abline(h=bestpast$theta[1],col="green");abline(h=true_theta[1],col="orange");abline(v=optloc,col="blue");  abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,2],type="l",ylab="2",     ylim=c( min(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])-1,max(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])+1)); abline(h=bestpast$theta[2],col="green");abline(h=true_theta[2],col="orange");abline(v=optloc,col="blue");  abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,3],type="l",ylab="sigma", ylim=c( min(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])-1,max(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])+1)); abline(h=bestpast$theta[3],col="green");abline(h=true_theta[3],col="orange");abline(v=optloc,col="blue");  abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,4],type="l",ylab="cap",   ylim=c( min(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])-1,max(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])+1)); abline(h=bestpast$theta[4],col="green");abline(h=true_theta[4],col="orange");abline(v=optloc,col="blue");  abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,5],type="l",ylab="fit",   ylim=c(0, max(true_theta[5]*2,5*min(ttrace[,5],na.rm=TRUE))));abline(h=true_theta[5],col="green");    abline(h=min(ttrace[,ncol(trace)],na.rm=TRUE),col="orange");     abline(v=optloc,col="blue");
    
    
    dev.flush(level = 1L)
    par(mfrow=c(1,1))
    },
    error = function(err) {
      print(paste("plotting error:  ",err))
    }
    ,finally={})
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

#f4<-function(x,prec=Inf,noiseseed){  1.5*(x[1]-0.987654321)*(x[2]-0.123456789)+(x[1]-0.987654321)^2+(x[2]-0.123456789)^2+(x[4]-0.123456789)^2 + x[5]^2+ (x[3]-0.123456789)^2+rnorm(1)/(prec)}
# f2<-function(x,prec=Inf){  1.5*(x[1]-0.987654321)*(x[2]-0.123456789)+(x[1]-0.987654321)^2+(x[2]-0.123456789)^2+rnorm(1)/(prec)}
# f1<-function(x,prec=10){  x^2+rnorm(1)/(prec)}
#iter=1000
#seed=round(runif(1)*1000)
#set.seed(seed)

#aaa<-random_walk_cooling(f4,c(1,1,1,1,1),iter=1000,minstep = 0.01, precschedule=function(iter,maxiter){1+{maxiter-iter}^2},upper=10,lower=-10,
#                                                 true_theta = c(0.987654321,rep(0.123456789,3),0))
#print(aaa$theta)
#x<-(c(rowSums(abs(sweep(aaa$trace,2,aaa$trace[nrow(iter),]) )),0))
#plot(x[20:iter])
