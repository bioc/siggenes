denspr<-function(x,n.interval=NULL,df=7){
	require(splines)
	if(is.null(n.interval))
		n.interval<-floor(sqrt(length(x)))
	breaks<-seq(min(x),max(x),length=n.interval+1)
	valHist<-hist(x,breaks=breaks,plot=FALSE)
	center<-valHist$mids
	counts<-valHist$counts
	ids<-which(counts>0)
	center<-center[ids]
	tmp<-ns.out<-ns(center,df=df)
	class(tmp)<-"matrix"
	mat<-data.frame(Number=counts[ids],tmp)
	glm.out<-glm(Number~.,data=mat,family=poisson)
	scale<-sum(diff(breaks)*counts)
	newx<-predict(ns.out,x)
	class(newx)<-"matrix"
	preds<-predict(glm.out,data.frame(newx),type="response")
	preds/scale
}

	


