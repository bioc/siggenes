pairt.samp.transform<-function(mat){
	if(!all(mat%in%c(1,-1)))
		stop("Some of the values in mat.samp are neither 1 nor -1.")
	vec.samp<-as.vector(t(mat))
	n.samp<-length(vec.samp)
	vec.new<-numeric(2*n.samp)
	vec.new[2*(1:n.samp)-1]<-vec.samp>0
	vec.new[2*(1:n.samp)]<-vec.samp<0
	mat.samp<-matrix(vec.new,nrow(mat),2*ncol(mat),byrow=TRUE)
	mat.samp
}
	 
 

