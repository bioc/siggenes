formatSAM<-function(x,digits=3){
	x<-format.pval(x,digits=digits,eps=0)
	x[as.numeric(x)==0]<-0
	x
}

pretty.mat.fdr<-function(x,digits=3){
	x<-as.data.frame(x)
	x[,"FDR"]<-formatSAM(x[,"FDR"],digits=digits)
	if(any(colnames(x)=="False")){
		x[,"False"]<-sapply(x[,"False"], 
			function(x) if(x>10^-digits) round(x,digits) else format(x,digits=3))
		idnr<-c(1,3,5)
	}
	else
		idnr<-c(1,3)
	x[,-idnr]<-round(x[,-idnr],digits)
	x
}

pretty.mat.sig<-function(x,digits=3){
	if(any(colnames(x)=="rawp")){
		x[,"rawp"]<-formatSAM(x[,"rawp"],digits=digits)
		x[,"q.value"]<-formatSAM(x[,"q.value"],digits=digits)
	}
	else
		x[,"local.fdr"]<-formatSAM(x[,"local.fdr"],digits=digits)
	format(x,digits=digits)
}