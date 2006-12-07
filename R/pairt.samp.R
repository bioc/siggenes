pairt.samp<-function(n.pair){
	mat<-matrix(0,2^n.pair,n.pair)
	for(i in 1:n.pair)
		mat[,i]<-rep(rep(c(1,-1),e=2^(n.pair-i)),2^(i-1))
	mat
} 
 

