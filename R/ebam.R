`ebam` <- function(x, cl, method=z.ebam, delta=.9, which.a0=NULL, control=ebamControl(),
		gene.names=dimnames(x)[[1]], ...){
	xclass<-class(x)
	if(!xclass%in%c("FindA0","ExpressionSet","matrix","data.frame", "list"))
		stop("x must be an object of class FindA0, ExpressionSet,\n",
			"matrix, data.frame, or list.")
	FUN<-match.fun(method)
	if(missing(cl)){
		if(xclass=="FindA0"){
			out<-ebamA0(x,which.a0=which.a0)
			chip.name<-x@chip
		}
		else{
			chip.name <- ""
			usedFun <- deparse(substitute(method))
			if(usedFun %in% c("chisq.ebam", "trend.ebam")){
				if(xclass!="list")
					stop("x must be a list.")
				out <- FUN(x, ...)
			}
			else
				stop("cl needs to be specified.")
		}
	}
	else{
		if(is(x,"ExpressionSet")){
			# requireNamespace("affy", quietly=TRUE)
			chip.name <- x@annotation
			if(is.character(cl) & length(cl)<=2)
				cl<-pData(x)[,cl]
			x<-exprs(x)
		}
		else
			chip.name<-""
		if(is.factor(cl))
			cl<-as.character(cl)
		if(ncol(x)!=length(cl))
			stop("The number of columns of data must be equal to the length of cl.")
		out<-FUN(x,cl,...)
	}
	check.out<-checkFUNout(out)
	if(is.na(control$p0))
		p0 <- pi0.est3(out, control$p0.estimation, exact=check.out$exact,
			lambda=control$lambda, ncs.value=control$ncs.value,
			use.weights=control$use.weights)
	else
		p0 <- control$p0
	if(!is.null(gene.names)){
		names(out$z)<-names(out$ratio)<-substring(gene.names,1,50)
		if(length(check.out$vec.pos)!=0)
			names(check.out$vec.pos)<-names(check.out$vec.neg)<-names(out$z)
	}
	posterior<-1-p0*out$ratio
	posterior[posterior<0]<-0
	mat<-compNumber(out$z,posterior,p0,check.out$B,delta=delta,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg)	
	new("EBAM",z=out$z,posterior=posterior,p0=p0,local=1-posterior,mat.fdr=mat,
		a0=check.out$a0,mat.samp=check.out$mat.samp,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg,msg=check.out$msg,chip=chip.name)		
		
}

