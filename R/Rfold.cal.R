Rfold.cal<-function(mat,cl,unlog=TRUE,R.fold=1,use.dm=FALSE){
	if(R.fold<=0)
		stop("R.fold must be larger than 0.")
	if(R.fold<1)
		R.fold<-1/R.fold
	if(!use.dm && unlog)
		mat<-2^mat
	mean0<-rowMeans(mat[,cl==0])
	mean1<-rowMeans(mat[,cl==1])
	if(!use.dm){
		mean0[mean0<=0]<-NA
		mean1[mean1<=0]<-NA
	}
	fold<-if(!use.dm) mean1/mean0 else 2^(mean1-mean0)
	fulfill<-numeric(length(fold))
	fulfill[fold>=R.fold | fold<=1/R.fold]<-1
	cbind(fold=fold,fulfill=fulfill)
}
	 
 

