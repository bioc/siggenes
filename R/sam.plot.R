# Copyright (C) 2003 Holger Schwender

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


sam.plot<-function(sam.out,delta,q.values=TRUE,R.fold=TRUE,R.unlog=TRUE,
		na.rm=FALSE,file.out=NA,gene.names=NULL){
	use.numbers<-sam.out$use.numbers
	rand<-sam.out$rand
	data<-eval(sam.out$X.name)
	if(R.unlog)
		data<-2^data
	if(!is.null(gene.names))
		data<-cbind(data,gene.names)
	cgn<-ifelse(is.null(gene.names),NA,ncol(data))
	tmp<-sam.plot.old(sam.out,delta,data,q.values=q.values,R.fold=R.fold,na.rm=na.rm,
		file.out=file.out,col.gene.name=cgn,use.numbers=use.numbers,rand=rand)
	invisible(structure(list(vec.fdr=tmp$vec.fdr,sam.output=tmp$sam.output,
		row.sig.genes=tmp$row.sig.genes)))
}
	