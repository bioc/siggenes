cat.ebam<-function(x,y,approx=FALSE,B=100,n.split=1,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,fast=FALSE,n.interval=NULL,df.ratio=3,
		df.dens=NULL,knots.mode=NULL,type.nclass="wand",rand=NA){
	x<-as.matrix(x)
	if(any(is.na(x)))
		stop("No NAs allowed.")
	if(check.for.NN){
		if(any(x=="NN"))
			stop("No NNs allowed.")
	}
	if(ncol(x)!=length(y))
		stop("The number of columns of x must be equal to the length of y.")
	if(!is.null(lev))
		x<-recodeLevel(x,lev)
	if(mode(x)!="numeric")
		mode(x)<-"numeric"
	if(any(is.na(x)))
		stop("x must contain integers.")
	n.cat<-max(x)
	if(any(!x%in%1:n.cat))
		stop("x must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%x))
		stop("Some of the values between 1 and ",n.cat," are not in x.")
	if(n.split<=0)
		stop("n.split must be at least 1.")
	if(n.split==1)
		stats<-chisqClass(x,y,n.cat)    #,check=check.levels)
	else{
		if(!approx)
			stop("Currently, splitting the variables is not supported in the permutation case.")
		stats<-chisqClassSplitted(x,y,n.cat,n.split)    #,check=check.levels)
	}
	n.cl<-length(unique(y))
	msg<-c("EBAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	if(approx){
		fail.out<-cat.null.approx2(stats,n.cl,n.cat,n.interval=n.interval,df.dens=df.dens,
			knots.mode=knots.mode,type.nclass=type.nclass)
		return(list(z=stats,ratio=fail.out$ratio,vec.pos=fail.out$vec.pos,
			vec.neg=fail.out$vec.neg,msg=msg))
	}
	if(is.null(n.interval))
		n.interval<-139
	out<-getSuccesses(stats,n.interval=n.interval)
	mat.samp<-setup.mat.samp(y,"f",B=B,B.more=B.more,B.max=B.max,rand=rand)
	fail.out<-compFailure(x,mat.samp,stats,out$interval,n.subset=n.subset,
		fast=fast,n.cat=n.cat)
	ratio<-compRatio(out$center,out$success,fail.out$vec.fail,df=df.ratio,z=stats)$ratio
	if(fast)
		return(list(z=stats,ratio=ratio,success=out$success,failure=fail.out$vec.fail,
			center=out$center,mat.samp=mat.samp,msg=msg))
	else
		return(list(z=stats,ratio=ratio,vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B,mat.samp=mat.samp,msg=msg))
}

