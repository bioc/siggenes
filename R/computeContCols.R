"computeContCols" <-
function(x,CL,vec.ncl,n.obs,pam=FALSE){
	x<-x%*%CL
	r<-rowSums(x)
	#if(any(r==0))
	#	stop("All rows of data must have the same number of categories.")
	tmp<-r%*%t(vec.ncl)/n.obs
	#rowSums((x-tmp)^2/tmp)
	if(pam)
		return(list(N.obs=x,N.exp=tmp,mat.stat=x*x/tmp))
	rowSums(x*x/tmp)
}

