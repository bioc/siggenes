`checkInitialDelta` <-
function(object,initial=NULL,isSAM=TRUE,prec=6){
	if(isSAM)
		dists<-abs(sort(object@d)-object@d.bar)
	if(!is.null(initial)){
		if(length(initial)!=2 | !is.numeric(initial))
			stop("initial must be a numeric vector of length 2.",call.=FALSE)
		if(initial[1]>=initial[2])
			stop("The first value of initial must be smaller than the second value.",call.=FALSE)
		if(any(initial<=0))
			stop("Both values of initial must be larger than 0.",call.=FALSE) 
		if(!isSAM && any(initial>1))
			stop("For EBAM, all values of initial must be smaller than or equal to 1.",call.=FALSE)
		if(isSAM){
			if(initial[1]>=max(dists))
				stop("The first value of initial must be smaller than\n","the largest ",
					"value of Delta,","i.e. the maximum absolute distance between\n",
					"the observed and the corresponding expected test statistics.",
					call.=FALSE)
		}
		return(round(seq(initial[1],initial[2],le=10),prec))
	}
	if(isSAM){
		tmp<-range(dists)
		return(round(seq(max(0.1,tmp[1]),max(1,tmp[2]),le=10),prec))
	}
	tmp<-max(0.1,min(object@posterior))
	round(seq(tmp,max(object@posterior),le=10),prec)
}

