

bootstrap_parameter_estimate<-function(par,bootstrapseed,optimizer,drawcommand,drawcommand_control=list(),vdata,...,outcome="y") {
  
  datahash<-digest::digest(vdata)
  optimizerhash<-digest::digest(deparse(optimizer))
  
  dir.create("estimates")
  dir.create(paste0("estimates/",optimizerhash,sep=""))
  dir.create(paste0("estimates/",optimizerhash,"/",datahash,sep=""))
  
  dir.create("estimates")
  dir.create(paste0("estimates/",datahash,sep=""))
  dir.create(paste0("estimates/",datahash,"/",bootstrapseed,sep=""))
  dir.create(paste0("estimates/",datahash,"/",bootstrapseed,"/",optimizerhash,sep=""))

  storedestimates_file<-paste0("estimates/",optimizerhash,"/",datahash,"/boostrapestimates.R",sep="")
  storedestimates_file_new<-paste0("estimates/",datahash,"/",bootstrapseed,"/",optimizerhash,"/bootstrapestimate.Rda",sep="")
  
  callkey<-paste(sep = "",bootstrapseed,digest::digest(c(par,digest::digest(list(...)))))
  
  #open file if exists, use empty list otherwise
  boostrapestimates<-list()
  #if it exists in new files
  tryCatch(boostrapestimates<-readRDS(storedestimates_file_new), error = function(e) {})
  if (callkey %in% names(boostrapestimates))
  {
    thisest<-boostrapestimates[[callkey]]
    if ("val" %in% names(thisest)) {
      if (thisest$val>100) {
        cat("-")
        boostrapestimates[[callkey]]<-NULL
        saveRDS(boostrapestimates,storedestimates_file_new)
      } else if (thisest$val>10) {
        cat(thisest$val)
        cat("*")
      } else {
        cat("*")
      }
    }
    return(thisest)
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
  saveRDS(boostrapestimates,storedestimates_file_new)
  return(bootest)
}
