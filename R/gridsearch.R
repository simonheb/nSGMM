colMadss<-function(x) {colMeans(abs(sweep(x,2,colMeans(x))))}
#library(scales)
#library(mlrMBO)
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
                              #iter=1000,#how many iterations? =>rounds-currentround
                              #maxiter=0, #where did we start? =>rounds 
                              #tempresetiter=NULL, #when did we reset temperature? => rounds-lastreheat
                              rounds=1000, #
                              currentround=0, #
                              lastreheat=0,
                              momentum=0,#momentum direction (usually not specified)
                              momentumpersistence=0.5,#speed at which momentum decays
                              trace=NULL, #initial return value
                              stepsize_trace=NULL, 
                              stepsize_blowup=2,
                              stepsize_decay=2,
                              stepsize_rekindling=1,
                              minstep=0.01, #mimum stepsize
                              upper,  lower, #limits (important!)
                              precschedule,#a function that returns a integer that is passed as prec to fun
                              true_theta=NULL,#for plotting
                              noiseseed=1,
                              roundstopartial=200,
                              superverbose=FALSE,
                              current=NULL,currentvalueage=10,
                              lastplottet=Sys.time()-60,
                              bestpast=NULL,
                              true_g=NULL#the value a restet will use
) {
  if(rounds==currentround) {
    
    cat("\n")
    f1<-fun(bestpast$theta,prec=4*precschedule(rounds,rounds),noiseseed=noiseseed,...)
    f2<-fun(theta,prec=4*precschedule(rounds,rounds),noiseseed=noiseseed,...)
    if (f2>f1) {
      cat("returning something I found along the way",bestpast$currentround,"\n")
      #print(bestpast$theta)
      return(list(theta=bestpast$theta,
                  trace=trace,
                  value=f1,
                  bestpast=bestpast,
                  true_g=true_g))
    }
    #print(theta)
    return(list(theta=theta,
                trace=tail(trace,bestpast$currentround),
                value=f2,
                bestpast=bestpast,
                true_g=true_g))
  }
  
  theta<-pmin(theta,upper)
  theta<-pmax(theta,lower)
  
  
  
  if (is.null(stepsize)) stepsize <-(upper-lower)/100
  if (!any(is.null(true_theta)) & is.null(true_g)) {
    true_g<- fun(th=true_theta,prec=4*precschedule(1,rounds),noiseseed=noiseseed,...) 
  }  else {
    true_g<- -1
  } 
  
  if(currentround%%roundstopartial==50) { #every 100 rounds solve for the partial min
    theta<-partialoptim(fun,theta,prec=precschedule(rounds,rounds),lower=lower,upper=upper,steps=4,maxrounds=precschedule(rounds,rounds)/4,noiseseed=noiseseed,...)
    currentvalueage<-Inf #make sure this gets re-computed.
  }
  #decay momentum
  momentum<-momentum*momentumpersistence
  
  step<-rnorm(length(theta))*stepsize+momentum
  newtheta<-theta+step
  theta<-pmax(pmin(theta,upper),lower)
  newtheta<-pmax(pmin(newtheta,upper),lower)
  
  if(currentvalueage>=40 | is.null(current)) { #after recycling the same value 10 times, let's  use a new draw so that we are not getting stuck in a lucky draw from g
    old<-max(0.000001,fun(theta,prec=precschedule(rounds,rounds),noiseseed=noiseseed,...))
    currentvalueage<-0
    cat("=")
  }else {
    old<-current
  }
  #cat(old,"=",current,"\n")
  
  if (old==Inf & 0==currentround) {browser();stop("starting at a point that evaluates to Inf is unlikely to yield a good result")}
  new<-max(0.000001,fun(newtheta,prec=precschedule(rounds,rounds),...,noiseseed=noiseseed))
  if(old==Inf & new==Inf) {old<-new<-999999999} #because sometimes both are Inf
  temp<-runif(1)
  if(temp<(min(old/new,1)^((-lastreheat+currentround)/warmuplength))) {
    ##sometimes this runs into oblivion again, whenever this happens, get it back home
    #if (iter<600 & new>60*max(bestpast$fit,new) | iter<300 & new>20*max(bestpast$fit,new) | new>200*max(bestpast$fit,new)) { #can only be true if bestpas!=NULL
    #  newtheta<-bestpast$theta
    #  new<-bestpast$fit
    #  cat("*")
    #}
    effectivestep<-newtheta-theta
    theta<-newtheta
    current<-new
    currentvalueage<-0
    stepsize<-stepsize_blowup * stepsize + #blow up stepsize to be more probing, 
      stepsize_rekindling*1/(currentround+1) #add some extra searching early on
    if (currentround>2) stepsize<-pmin(stepsize,(upper-lower)/5) #prevent explosion, e.g. because prolonged warmup
    #only build up momentum if the step as an improvement
    #if (new<old) 
    momentum<-momentum+effectivestep
    cat("|")
    if (new<old | is.null(bestpast))  { #if the new value is better than old
      if (!is.null(bestpast))
        bestpast_fit<-bestpast$fit
      else 
        bestpast_fit<-Inf
      if (new<bestpast_fit) {
        bestpast<-list(theta=newtheta,currentround=currentround,fit=new)
      } else { #hope is, this would happen rarely  #better than old, but worse than the pastbest, maybe the pastbest was too optimist, so let's re-evaluate it so that we won't be stuck with a lucky draw
        cat("o")
        bestpast$fit<-max(0.000001,mean(c(bestpast$fit,fun(bestpast$theta,prec=precschedule(currentround,rounds),noiseseed=noiseseed,...))))
      }
    }
  } else {
    stepsize<-pmax(stepsize/stepsize_decay,minstep)
    cat(".")
    current<-old
    currentvalueage<-currentvalueage+1
  }
  
  trace=rbind(trace,c(theta,current))
  stepsize_trace=rbind(stepsize_trace,c(stepsize))
  if((rounds-currentround)%%reheatineval==0) {
    #reset the chain
    if (!is.null(bestpast)) {
      if (bestpast$currentround<currentround & rounds-50  > currentround) {
        bestpast$fit<-max(0.000001,fun(bestpast$theta,prec=precschedule(currentround,rounds),noiseseed=noiseseed,...)) #first, update
        if (bestpast$fit<current) { #if the past is better, reset.
          cat("[",bestpast$currentround,"]") 
          lastreheat<-currentround
          theta<-bestpast$theta
          stepsize<-head(colSd(trace[ceiling(0.9*nrow(trace)):nrow(trace),]),length(theta))
        }
      }
    }
  }
  
  
  if((rounds-currentround-1)%%500==0|rounds<=currentround+1) {
    if (is.null(true_theta)) {
      plot_trace(trace,c(theta,true_g),bestpast,rounds,stepsize_trace=stepsize_trace)
    } else {
      plot_trace(trace,c(true_theta,true_g),bestpast,rounds,stepsize_trace=stepsize_trace)
    } 
    lastplottet<-Sys.time()
    cat("\n[cround",currentround,":c(",paste0(round(theta,3),collapse=","),")=",current,"; pastbest@",bestpast$currentround,", was",bestpast$fit,"with theta=c(",paste0(round(bestpast$theta,3),collapse=","),")]\n")
    
  }
  return(random_walk_cooling(fun,theta,...,stepsize=stepsize,stepsize_rekindling=stepsize_rekindling,noiseseed=noiseseed,stepsize_decay=stepsize_decay,stepsize_blowup=stepsize_blowup,stepsize_trace=stepsize_trace,currentround=currentround+1,lastreheat=lastreheat,momentum=momentum,minstep=minstep,momentumpersistence=momentumpersistence,trace=trace,rounds=rounds,reheatineval=reheatineval,warmuplength=warmuplength,upper=upper,lower=lower,true_theta=true_theta,precschedule=precschedule,bestpast=bestpast,current=current,currentvalueage=currentvalueage,lastplottet=lastplottet,true_g=true_g))
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
    if(lower[par]==upper[par]) {
      cat("c")
    } else {
      xs_lower<- c( (1-((1:steps)/steps)^2)*(theta[par]-lower[par])+lower[par])
      xs_upper<- c( ((1:steps)/steps)^2*(upper[par]-theta[par])+theta[par])
      xs<-unique(c(xs_lower,xs_upper))
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
      } else {
        cat(">")
      }
    }
    par_i<-par_i+1
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
add.error.bars <- function(X,Y,SE,w,col=1){
  X0 = X; Y0 = (Y-SE); X1 =X; Y1 = (Y+SE);
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col);
}
plot_trace<-function(trace,true_theta,bestpast,rounds,stepsize_trace) {
  tryCatch({
    par(mfrow=c(ncol(trace),1))
    hideshare<-0.1
    ttrace<-trace
    if (nrow(trace)>100) 
      ttrace[1:ceiling(hideshare*nrow(trace)), 1:(ncol(trace))]<-NA
    
    optloc<-bestpast$currentround#which.min2(ttrace[,ncol(trace)]);
    ttrace[ttrace[,ncol(trace)]>max(true_theta[5]*2,7*min(ttrace[,ncol(trace)],na.rm=TRUE)),ncol(trace)]<-NA
    plot(ttrace[,1],type="l",ylab="1",     ylim=c( min(0,true_theta[1],tail(ttrace[,1],1),trace[optloc,1])-1,max(0,true_theta[1],tail(ttrace[,1],1),trace[optloc,1])+1));
      add.error.bars(1:nrow(stepsize_trace),ttrace[,1],stepsize_trace[,1],0,"cadetblue2");
      abline(h=bestpast$theta[1],col="green");abline(h=true_theta[1],col="orange");abline(v=optloc,col="blue");
      abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,2],type="l",ylab="2",     ylim=c( min(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])-1,max(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])+1));
      add.error.bars(1:nrow(stepsize_trace),ttrace[,2],stepsize_trace[,2],0,"cadetblue2");
      lines(ttrace[,2],type="l",ylab="2",     ylim=c( min(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])-1,max(0,true_theta[2],tail(ttrace[,2],1),trace[optloc,2])+1));
      abline(h=bestpast$theta[2],col="green");abline(h=true_theta[2],col="orange");abline(v=optloc,col="blue");
      abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,3],type="l",ylab="sigma", ylim=c( min(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])-1,max(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])+1));
      add.error.bars(1:nrow(stepsize_trace),ttrace[,3],stepsize_trace[,3],0,"cadetblue2");
      lines(ttrace[,3],type="l",ylab="sigma", ylim=c( min(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])-1,max(0,true_theta[3],tail(ttrace[,3],1),trace[optloc,3])+1));
      abline(h=bestpast$theta[3],col="green");abline(h=true_theta[3],col="orange");abline(v=optloc,col="blue");
      abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,4],type="l",ylab="cap",   ylim=c( min(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])-1,max(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])+1));
      add.error.bars(1:nrow(stepsize_trace),ttrace[,4],stepsize_trace[,4],0,"cadetblue2");
      lines(ttrace[,4],type="l",ylab="cap",   ylim=c( min(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])-1,max(0,true_theta[4],tail(ttrace[,4],1),trace[optloc,4])+1));
      abline(h=bestpast$theta[4],col="green");abline(h=true_theta[4],col="orange");abline(v=optloc,col="blue");
      abline(h = 0, lty = 2,col="gray")
    plot(ttrace[,5],type="l",ylab="fit",   ylim=c(0, max(true_theta[5]*2,5*min(ttrace[,5],na.rm=TRUE))));
      abline(h=true_theta[5],col="green");
      abline(h=min(ttrace[,ncol(trace)],na.rm=TRUE),col="orange");
      abline(v=optloc,col="blue");

    
    dev.flush(level = 1L)
    par(mfrow=c(1,1))
  },
  error = function(err) {
    print(paste("plotting error:  ",err))
  }
  ,finally={})
}
