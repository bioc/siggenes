# Copyright (C) 2003 Holger Schwender

ebam.wilc<-function(data,cl,delta=.9,p0=NA,ties.rand=TRUE,zero.rand=TRUE,gene.names=NULL,
		R.fold=TRUE,R.unlog=TRUE,file.out=NA,na.rm=FALSE,rand=NA){
	xy.out<-xy.cal(cl,T,T)
	x<-xy.out$x
	y<-xy.out$y
	paired<-xy.out$paired
	data2<-if(R.unlog) 2^data else data
	if(!is.null(gene.names))
		data<-cbind(data,gene.names)
	cgn<-ifelse(is.null(gene.names),NA,ncol(data))
	tmp<-ebam.wilc.old(data,x,y,paired=paired,delta=delta,p0=p0,ties.rand=ties.rand,zero.rand=zero.rand,
		col.gene.name=cgn,R.fold=R.fold,R.dataset=data2,file.out=file.out,rand=rand,
		na.rm=na.rm)
	invisible(structure(list(nsig=tmp$nsig,false=tmp$false,fdr=tmp$fdr,ebam.out=tmp$ebam.out,
		mat.out=tmp$mat.out,p0=tmp$p0,glm.out=tmp$glm.out,f.x=tmp$f.x,f.null=tmp$f.null,
		y.wilc=tmp$y.wilk,ebam.output=tmp$ebam.output,row.sig.genes=tmp$row.sig.genes)))
}