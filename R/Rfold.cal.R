Rfold.cal<-function(mat,cl,unlog=TRUE,R.fold=1){
	if(R.fold<=0)
		stop("R.fold must be larger than 0.")
	if(R.fold<1)
		R.fold<-1/R.fold
	if(unlog)
		mat<-2^mat
	uni.cl<-sort(unique(cl))
	g1<-if(length(uni.cl)==2) which(cl==uni.cl[1]) else which(cl<0)
	mean.g1<-rowMeans(mat[,g1])
	mean.g1[mean.g1<=0]<-NA
	mean.g2<-rowMeans(mat[,-g1])
	mean.g2[mean.g2<=0]<-NA
	fold<-mean.g2/mean.g1
	fulfill<-numeric(length(fold))
	fulfill[fold>=R.fold | fold<=1/R.fold]<-1
	cbind(fold=fold,fulfill=fulfill)
}
	 
 

