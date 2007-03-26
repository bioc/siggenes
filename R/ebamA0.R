`ebamA0` <-
function(x,which.a0=NULL){
	if(is.null(which.a0)){
		a0<-x@suggested
		which.a0<-which(x@vec.a0==a0)
	}
	else{
		if(!which.a0%in%(1:length(x@vec.a0)))
			stop("which.a0 must be an integer between 1 and ",length(x@vec.a0),".")
		a0<-x@vec.a0[which.a0]
	}
	ratio<-compRatio(x@mat.center[,which.a0],x@mat.success[,which.a0],x@mat.failure[,which.a0],
		df=x@df.ratio,z=x@mat.z[,which.a0])$ratio
	list(z=x@mat.z[,which.a0],ratio=ratio,a0=a0,success=x@mat.success[,which.a0],
		failure=x@mat.failure[,which.a0],center=x@mat.center[,which.a0],
		mat.samp=x@mat.samp,msg=x@msg)
}

