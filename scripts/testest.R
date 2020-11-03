estimates1<-tryCatch(readRDS("estimates-m4m6m7-test-convergence-optima22"),error=function(e){NULL},finally={})
estimates<-tryCatch(readRDS("estimates-m4m6m7-test-convergence-optima-restart"),error=function(e){NULL},finally={})


estimates<-merge(estimates1,estimates2,by=c("village","type"))
estimates$impliedconstr.x<-1/estimates$theta.4.x
estimates$impliedconstr.y<-1/estimates$theta.4.y

for (i in 1:10) {

g(c(-3.3524870,  4.354150,  0.9690552, invkappatransformation(1/5.166298e+03)),kinship=vdata[[village]][["m8am8bm8c"]],
  income=vdata[[village]][["income"]]+1, keep=keep_with_stars,
  transfers=vdata[[village]][[transfernet]],
  distance=vdata[[village]][["distance"]],noiseseed=3998,prec=2000,vcv=var(mcpp),verbose=F)


g(c(-0.3203979,  2.068022,  0.1239922, invkappatransformation(1/2.915678e-01)),kinship=vdata[[village]][["m8am8bm8c"]],
  income=vdata[[village]][["income"]]+1, keep=keep_with_stars,
  transfers=vdata[[village]][[transfernet]],
  distance=vdata[[village]][["distance"]],noiseseed=3998,prec=2000,vcv=var(mcpp),verbose=F)
}