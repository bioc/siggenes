cat.ebam<-function(data,cl,approx=FALSE,B=100,n.split=1,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,fast=FALSE,n.interval=NULL,df.ratio=3,
		df.dens=NULL,knots.mode=NULL,type.nclass="wand",rand=NA){
	data<-as.matrix(data)
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.")
	if(check.for.NN && any(data=="NN" | data=="NoCall")){
		if(approx)
			data[data=="NN" | data=="NoCall"] <- NA
		else
			stop("No missing calls allowed.")
	}
	if(any(is.na(data))){
		if(approx)
			withNA <- TRUE
		else
			stop("No missing values allowed in data when approx = FALSE.")
	}
	else
		withNA <- FALSE
	if(!is.null(lev))
		data<-recodeLevel(data,lev)
	if(mode(data)!="numeric")
		mode(data)<-"numeric"
	n.cat<-max(data, na.rm=TRUE)
	if(any(!data%in%c(1:n.cat, NA)))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	if(n.split<=0)
		stop("n.split must be at least 1.")
	if(n.split==1)
		stats<-chisqClass(data,cl,n.cat,withNA=withNA)
	else{
		if(!approx)
			stop("Currently, splitting the variables is not supported in the permutation case.")
		stats<-chisqClassSplitted(data,cl,n.cat,n.split,withNA=withNA)
	}
	n.cl<-length(unique(cl))
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

