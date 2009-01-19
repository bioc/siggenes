rowWilcoxon <- function(X, y, rand=NA){
	cat("Currently under testing!","\n\n")
	if(any(!y%in%c(0,1)))
		stop("y must consist of zeros and ones.")
	if(all(y==0))
		stop("None of the class labels is equal to one.")
	if(any(y==0)){
		out<-rowRanksWilc(X, y)
		return(rowSums(out))
	}
	if(any(X==0)){
		if(!is.na(rand))
			set.seed(rand)
		warning("There are ",sum(X==0)," observations/pairs having a value/difference of zero.",
			"\n","These values/differences are randomly set to either 1e-06 or -1e-06.",
			call.=FALSE)
		X[X==0]<-sample(c(1e-06,-1e-06),sum(X==0),rep=TRUE)
	}
	out <- rowRanksWilc(abs(X), y)
	rowSums(out*(X>0))
}

