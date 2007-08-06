cat.stat<-function(data,cl,B=100,approx=FALSE,n.split=1,check.levels=TRUE,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,rand=NA){
	data<-as.matrix(data)
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(check.for.NN){
		if(any(data=="NN"))
			stop("No NNs allowed.")
	}
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(!is.null(lev))
		data<-recodeLevel(data,lev)
	if(mode(data)!="numeric")
		mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("data must contain integers.")
	n.cat<-max(data)
	if(any(!data%in%1:n.cat))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	if(n.split<=0)
		stop("n.split must be at least 1.")
	if(n.split==1)
		stats<-chisqClass(data,cl,n.cat,check=check.levels)
	else{
		if(!approx)
			stop("Currently, splitting the variables is not supported in the permutation case.")
		stats<-chisqClassSplitted(data,cl,n.cat,n.split,check=check.levels)
	}
	n.cl<-length(unique(cl))
	if(approx){
		null.out<-cat.null.approx(stats,n.cl,n.cat)
		mat.samp<-matrix(numeric(0))
	}
	else{
		mat.samp<-setup.mat.samp(cl,"f",B=B,B.more=B.more,B.max=B.max,
			rand=rand)
		null.out<-cat.null(data,mat.samp,stats,n.subset,n.cat)
	}
	msg<-c("SAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	structure(list(d=stats,d.bar=null.out$d.bar,p.value=null.out$p.value,
		vec.false=null.out$vec.false,discrete=FALSE,s=numeric(0),
		s0=numeric(0),mat.samp=mat.samp,msg=msg,fold=numeric(0)))
}

