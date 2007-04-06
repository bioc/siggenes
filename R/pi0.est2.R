`pi0.est2` <-
function(fail,success,z,type,lambda=NULL,ncs.value="max",use.weights=FALSE){
	if(!type%in%c("splines","interval"))
		stop("type must be either 'splines' or 'interval' in pi0.est2.")
	if(is.null(lambda)){
		lambda<-if(type=="splines") seq(0,0.95,0.05) else 0.5
	}
	w<-if(use.weights) 1-lambda else NULL
	z.order<-order(abs(z))
	z.fail<-cumsum(fail[z.order])
	z.success<-cumsum(success[z.order])
	n.fail<-sum(fail)
	B<-n.fail/sum(success)
	if(type=="interval"){
		num<-n.fail*(1-lambda)
		id<-which.min(abs(z.fail-num))
		pi0<-B*z.success[id]/(z.fail[id])
	}
	else{
		tmp<-(n.fail-z.fail[-length(z.fail)])/n.fail
		p<-rep(c(1,tmp),success[z.order])
		pi0<-pi0.est(p,lambda=lambda,ncs.value=ncs.value,ncs.weights=w)$p0
	}	
	pi0	
}

