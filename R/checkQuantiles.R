`checkQuantiles` <-
function(s,quan.a0=(0:10)/100,include.zero=TRUE){
	if(length(quan.a0)<2)
		stop("quan.a0 must be a vector containing at least two values.")
	if(any(quan.a0<0 | quan.a0>1))
		stop("The values in quan.a0 must be between 0 and 1.")
	if (any(round(100 * quan.a0, 10) != round(100 * quan.a0,0))){
        	warning("At least one alpha is not a percentile. Only the first two decimal digits", 
            		" are retained.")
        	quan.a0 <- signif(quan.a0, 2)
    	}
	quans<-quantile(s,quan.a0)
	if(include.zero)
		quans<-c(0,quans)
	quans
}

