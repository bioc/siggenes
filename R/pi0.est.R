pi0.est<-function(p,lambda=seq(0,.95,.05),ncs.value="max",ncs.weights=NULL){
	if(any(lambda>=1) || any(lambda<0))
		stop("lambda has to be between 0 and 1.")
	if(!ncs.value%in%c("paper","max"))
		stop("ncs.value has to be either paper oder max.")
	le.la<-length(lambda)
	ncs.pred<-ifelse(ncs.value=="paper",1,max(lambda))
	vec.p0<-numeric(le.la)
	for(i in 1:le.la)
		vec.p0[i]<-mean(p>lambda[i])/(1-lambda[i])
	if(le.la==1)
		return(list(p0=min(vec.p0,1),spline.out=NULL))
	spline.out<-smooth.spline(lambda,vec.p0,w=ncs.weights,df=3)
	p0<-min(predict(spline.out,ncs.pred)$y,1)
	if(p0<=0){
		warning("The spline based estimation of pi0 results in a non-positive value of pi0.",
			"\n","Therefore, pi0 is estimated by using lambda = 0.5.",call.=FALSE)
		p0<-vec.p0[lambda==.5]
		spline.out<-NULL
	}
	return(list(p0=p0,spline.out=spline.out))
} 
 

