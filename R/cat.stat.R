cat.statNew <- function(...){
	stop("\n","cat.stat has been removed. Please use one of the following methods instead.",
		"\n\n",
		"chisq.stat   Pearson's ChiSquare Statistic (was cat.stat)\n",
		"chisq.tab    Pearson's ChiSquare Statistic based on Tables\n",
		"catt.stat    Cochran-Armitage Trend Test\n",
		"catt.tab     Cochran-Armitage Trend Test based on Tables\n",
		"trend.stat   General Trend Test\n",
		"trend.tab    General Trend Test based on Tables\n",
		"\n", "Still to come:\n",
		"pccc.stat    Pearson's Corrected Contingency Coefficient\n",
		"pccc.tab     Pearson's Corrected Contingency Coefficient for Tables\n",
		"\n", "It is recommended to use *.tab if the number of samples is very large.\n")
}



cat.stat<-function(x,y,B=100,approx=FALSE,n.split=1,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,rand=NA){
	x<-as.matrix(x)
	if(any(is.na(x)))
		stop("No missing values allowed.")
	if(check.for.NN){
		if(any(x=="NN" | x=="NoCall"))
			stop("No missing calls allowed.")
	}
	if(ncol(x)!=length(y))
		stop("The number of columns of x must be equal to the length of y.")
	if(!is.null(lev))
		x<-recodeLevel(x,lev)
	if(mode(x)!="numeric")
		mode(x)<-"numeric"
	if(any(is.na(x)))
		stop("No missing values allowed in 'x'.")
	n.cat<-max(x)
	if(any(!x%in%1:n.cat))
		stop("x must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%x))
		stop("Some of the values between 1 and ",n.cat," are not in x.")
	if(n.split<=0)
		stop("n.split must be at least 1.")
	if(n.split==1)
		stats<-chisqClass(x,y,n.cat) #,check=check.levels)
	else{
		if(!approx)
			stop("Currently, splitting the variables is not supported in the permutation case.")
		stats<-chisqClassSplitted(x,y,n.cat,n.split) #,check=check.levels)
	}
	n.cl<-length(unique(y))
	if(approx){
		null.out<-cat.null.approx(stats,n.cl,n.cat)
		mat.samp<-matrix(numeric(0))
	}
	else{
		mat.samp<-setup.mat.samp(y,"f",B=B,B.more=B.more,B.max=B.max,
			rand=rand)
		null.out<-cat.null(x,mat.samp,stats,n.subset,n.cat)
	}
	msg<-c("SAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	structure(list(d=stats,d.bar=null.out$d.bar,p.value=null.out$p.value,
		vec.false=null.out$vec.false,discrete=FALSE,s=numeric(0),
		s0=numeric(0),mat.samp=mat.samp,msg=msg,fold=numeric(0)))
}

