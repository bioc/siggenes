`compRatio` <-
function(center,succ,fail,df=5,z=NULL){ 
	# requireNamespace("splines")
	n<-succ+fail
	ids<-which(n>0)
	if(is.null(z))
		z<-center
	tmp<-ns.out <- ns(center[ids], df)
	class(tmp)<-"matrix"
	mat<-data.frame(s=succ[ids],f=fail[ids],tmp)
	glm.out<-glm(cbind(s,f)~.,data=mat,family=binomial)
	newz<-predict(ns.out,z)
    	class(newz)<-"matrix"
	probs<-predict(glm.out,data.frame(newz),type="response")
	B<-sum(fail)/sum(succ)
	r<-(1-probs)/(B*probs)
	return(list(z=z,probs=probs,ratio=r))
}

