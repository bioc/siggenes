`pi0.est3` <-
function(fun.out,p0.estimation=c("splines","interval","adhoc"),exact=TRUE,lambda=NULL,
		ncs.value="max",use.weights=FALSE){
	type.p0<-match.arg(p0.estimation)
	if(type.p0=="adhoc")
		return(min(1/fun.out$ratio))
	if(!exact)
		return(pi0.est2(fun.out$failure,fun.out$success,fun.out$center,type.p0,
			lambda=lambda,ncs.value=ncs.value,use.weights=use.weights))
	if(is.null(lambda))
		lambda<-if(type.p0=="splines") seq(0,0.95,0.05) else 0.5
	w<-if(use.weights) 1-lambda else NULL
	p<-(fun.out$vec.pos+fun.out$vec.neg)/length(fun.out$vec.pos)
	pi0.est(p,lambda=lambda,ncs.value=ncs.value,ncs.weights=w)$p0
}

