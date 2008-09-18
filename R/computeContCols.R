"computeContCols" <-
function(x,CL,vec.ncl,n.obs){
	x<-x%*%CL
	r<-rowSums(x)
	tmp<-r%*%t(vec.ncl)/n.obs
	rowSums(x*x/tmp)
}

