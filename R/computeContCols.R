"computeContCols" <-
function(x,CL,vec.ncl,n.obs){
	x<-x%*%CL
	r<-rowSums(x)
	if(is.matrix(vec.ncl))
		tmp <- r * vec.ncl / n.obs
	else 
		tmp <- r %*% t(vec.ncl) / n.obs
	rowSums(x*x/tmp)
}

