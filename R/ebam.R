`ebam` <-
function(x,cl,method=z.ebam,delta=.9,which.a0=NULL,p0=NA,
		p0.estimation=c("splines","interval","adhoc"),lambda=NULL,ncs.value="max",
		use.weights=FALSE,gene.names=dimnames(x)[[1]],...){
	xclass<-class(x)
	if(!xclass%in%c("FindA0","ExpressionSet","matrix","data.frame"))
		stop("x must be an object of class FindA0, ExpressionSet, matrix, or data.frame.")
	if(xclass=="FindA0"){
		out<-ebamA0(x,which.a0=which.a0)
		chip.name<-x@chip
	}
	else{
		if(missing(cl))
			stop("cl must be specified if x is a matrix, a data frame, ",
				"or an ExpressionSet object.")
		if(is(x,"ExpressionSet")){
			require(affy,quietly=TRUE)
			chip.name<-annotation(x)
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
		FUN<-match.fun(method)
		out<-FUN(x,cl,...)
	}
	check.out<-checkFUNout(out)
	if(is.na(p0))
		p0<-pi0.est3(out,p0.estimation,exact=check.out$exact,lambda=lambda,
			ncs.value=ncs.value,use.weights=use.weights)
	if(!is.null(gene.names)){
		names(out$z)<-names(out$ratio)<-substring(gene.names,1,50)
		if(length(check.out$vec.pos)!=0)
			names(check.out$vec.pos)<-names(check.out$vec.neg)<-names(out$z)
	}
	posterior<-1-p0*out$ratio
	posterior[posterior<0]<-0
	mat<-compNumber(out$z,posterior,p0,check.out$B,delta=delta,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg)	
	new("EBAM",z=out$z,posterior=posterior,p0=p0,local=p0*out$ratio,mat.fdr=mat,
		a0=check.out$a0,mat.samp=check.out$mat.samp,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg,msg=check.out$msg,chip=chip.name)		
		
}

