chisq.stat <- function(data, cl, approx=NULL, B=100, n.split=1, check.for.NN=FALSE,
		lev=NULL, B.more=0.1, B.max=50000, n.subset=10, rand=NA){
	if(missing(cl))
		return(chisqStatMissing(data, approx=approx))
	if(is.data.frame(data))
		data <- as.matrix(data)
	if(!is.matrix(data))
		stop("data must be a matrix or a data.frame.")
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.")
	if(is.null(approx))
		approx <- FALSE
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
	if(any(!data %in% c(1:n.cat,NA)))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	if(n.split<=0)
		stop("n.split must be at least 1.")
	if(n.split==1)
		stats<-chisqClass(data,cl,n.cat, withNA=withNA) #,check=check.levels)
	else{
		if(!approx)
			stop("Currently, splitting the variables is not supported in the permutation case.")
		stats<-chisqClassSplitted(data,cl,n.cat,n.split, withNA=withNA) #,check=check.levels)
	}
	n.cl<-length(unique(cl))
	if(approx){
		null.out<-cat.null.approx(stats,n.cl,n.cat)
		mat.samp<-matrix(numeric(0))
	}
	else{
		mat.samp<-setup.mat.samp(cl-1,"f",B=B,B.more=B.more,B.max=B.max,
			rand=rand) + 1
		null.out<-cat.null(data,mat.samp,stats,n.subset,n.cat)
	}
	msg<-c("SAM Analysis for Categorical Data\n\n",
		paste("Null Distribution: ",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	structure(list(d=stats,d.bar=null.out$d.bar,p.value=null.out$p.value,
		vec.false=null.out$vec.false,discrete=FALSE,s=numeric(0),
		s0=numeric(0),mat.samp=mat.samp,msg=msg,fold=numeric(0)))
}

chisqStatMissing <- function(data, approx=NULL){
	require(scrime)
	if(!inherits(data, "list"))
		stop("data must be a list of matrices if cl is not specified.")
	if(is.null(approx))
		approx <- TRUE
	if(!approx)
		stop("Currently only approx=TRUE available if data is a list.")
	out <- if(length(data)==2) rowChisq2Class(data[[1]], data[[2]], sameNull=TRUE)
		else rowChisqMultiClass(listTables=data, sameNull=TRUE)
	n.var <- length(out$stats)
	d.bar <- qchisq(((1:n.var) - 0.5) / n.var, out$df)
	vec.false <- out$rawp * n.var
	msg <- c("SAM Analysis for Categorical Data\n\n",
		"Null Distribution: Approximation by ChiSquare Distribution\n\n")
	structure(list(d=out$stats, d.bar=d.bar, p.value=out$rawp, vec.false=vec.false,
		discrete=FALSE, s=numeric(0), s0=numeric(0), mat.samp=matrix(numeric(0)),
		msg=msg, fold=numeric(0)))
}
