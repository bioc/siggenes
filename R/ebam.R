# Copyright (C) 2003 Holger Schwender

ebam<-function(a0.out,a0=NA,p0=NA,delta=NA,local.bin=.1,gene.names=NULL,q.values=TRUE,
		R.fold=TRUE,R.unlog=TRUE,na.rm=FALSE,file.out=NA){
	data<-eval(a0.out$X.name)
	if(R.unlog) 
		data<-2^data
	if(!is.null(gene.names))
		data<-cbind(data,gene.names)
	cgn<-ifelse(is.null(gene.names),NA,ncol(data))
	tmp<-ebam.old(a0.out,data,a0=a0,p0=p0,delta=delta,local.bin=.1,col.gene.name=cgn,
		q.values=q.values,R.fold=R.fold,R.dataset=data,na.rm=na.rm,file.out=file.out)
	invisible(structure(list(mat.repeat=tmp$mat.repeat,optim.out=tmp$optim.out,
		mat.post.Z=tmp$mat.post.Z,ebam.out=tmp$ebam.out,FDR=tmp$FDR,a0=tmp$a0,
		mat.Z.unsorted=tmp$mat.Z.unsorted,row.sig.genes=tmp$row.sig.genes,p0=tmp$p0)))
}