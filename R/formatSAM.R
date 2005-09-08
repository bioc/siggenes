formatSAM<-function(x,digits=3){
	x<-format.pval(x,digits=digits,eps=0)
	x[as.numeric(x)==0]<-0
	x
}

pretty.mat.fdr<-function(x,digits=3){
	x<-as.data.frame(x)
	x[,"FDR"]<-formatSAM(x[,"FDR"],digits=digits)
	x[,"False"]<-sapply(x[,"False"], 
		function(x) if(x>10^-digits) round(x,digits) else format(x,digits=3))
	x[,-c(1,3,5)]<-round(x[,-c(1,3,5)],digits)
	x
}

pretty.mat.sig<-function(x,digits=3){
	x[,"rawp"]<-formatSAM(x[,"rawp"],digits=digits)
	x[,"q.value"]<-formatSAM(x[,"q.value"],digits=digits)
	format(x,digits=digits)
}