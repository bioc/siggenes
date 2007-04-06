cat.ebam<-function(data,cl,approx=FALSE,B=100,check.levels=TRUE,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,fast=FALSE,n.interval=139,df.ratio=3,
		df.glm=NULL,rand=NA){
	data<-as.matrix(data)
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(check.for.NN){
		if(any(data=="NN"))
			stop("No NNs allowed.")
	}
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(!is.null(lev)){
		for(i in 1:length(lev))
			data[data==lev[i]]<-i
	}
	if(mode(data)!="numeric")
		mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("data must contain integers.")
	n.cat<-max(data)
	if(any(!data%in%1:n.cat))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	stats<-chisqClass(data,cl,n.cat,check=check.levels)
	if(is.null(n.interval))
		n.interval<-ifelse(approx,floor(sqrt(length(stats))),139)
	n.cl<-length(unique(cl))
	msg<-c("EBAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	if(approx){
		fail.out<-cat.null.approx2(stats,n.cl,n.cat,n.interval=n.interval,df.glm=df.glm)
		return(list(z=stats,ratio=fail.out$ratio,vec.pos=fail.out$vec.pos,
			vec.neg=fail.out$vec.neg,msg=msg))
	}
	out<-getSuccesses(stats,n.interval=n.interval)
	mat.samp<-setup.mat.samp(cl,"f",B=B,B.more=B.more,B.max=B.max,rand=rand)
	fail.out<-compFailure(data,mat.samp,stats,out$interval,n.subset=n.subset,
		fast=fast,n.cat=n.cat)
	ratio<-compRatio(out$center,out$success,fail.out$vec.fail,df=df.ratio,z=stats)$ratio
	if(fast)
		return(list(z=stats,ratio=ratio,success=out$success,failure=fail.out$vec.fail,
			center=out$center,mat.samp=mat.samp,msg=msg))
	else
		return(list(z=stats,ratio=ratio,vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B,mat.samp=mat.samp,msg=msg))
}
