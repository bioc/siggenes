`checkA0` <-
function(s,a0=NULL,quan.a0=NULL){
	if(is.null(a0) & is.null(quan.a0))
		stop("Either a0 or quan.a0 must be specified. See the help of z.ebam.")
	if(!is.null(a0) & !is.null(quan.a0))
		stop("Only one of a0 or quan.a0 should be specified.")
	if(!is.null(a0)){
		if(length(a0)>1 | !is.numeric(a0))
			stop("a0 should be a numeric value (and not a vector of values).")
		if(a0<0)
			stop("a0 must be larger than or equal to 0.")
		return(a0)
	}
	if(length(quan.a0)>1 | !is.numeric(quan.a0))
		stop("quan.a0 must be a numeric value (and not, e.g., a vector).")
	if(quan.a0<0 | quan.a0>1)
		stop("quan.a0 must be between 0 and 1.")
	quantile(s,quan.a0)
}

