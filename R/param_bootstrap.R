

bootstrap_parameter_estimate<-function(par,bootstrapseed,optimizer,drawcommand,drawcommand_control=list(),vdata,...,outcome="y") {
  
  datahash<-digest::digest(vdata)
  optimizerhash<-digest::digest(deparse(optimizer))
  
  dir.create("estimates")
  dir.create(paste0("estimates/",optimizerhash,sep=""))
  dir.create(paste0("estimates/",optimizerhash,"/",datahash,sep=""))

  storedestimates_file<-paste0("estimates/",optimizerhash,"/",datahash,"/boostrapestimates.R",sep="")
  
  callkey<-paste(sep = "",bootstrapseed,digest::digest(c(par,digest::digest(list(...)))))
  
  #open file if exists, use empty list otherwise
  boostrapestimates<-list()
  tryCatch(boostrapestimates<-readRDS(storedestimates_file), error = function(e) {})
  
  if (callkey %in% names(boostrapestimates))
  {
    cat("*")
    #print(boostrapestimates[[callkey]])
    return(boostrapestimates[[callkey]])
  }
  ptm<-Sys.time()
  set.seed(bootstrapseed)
  #generate fake data
  vdata[[outcome]]<-
    do.call(drawcommand,c(list(par,vdata),drawcommand_control)) 
  #obtain estimates
  bootest<-optimizer(vdata,...)
  boostrapestimates[[callkey]]<-bootest
  boostrapestimates[[callkey]]$time<-round(as.numeric(Sys.time() - ptm,units="mins"))
  saveRDS(boostrapestimates,storedestimates_file)
  return(bootest)
}
