# Copyright (C) 2003 Holger Schwender

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


sam<-function(data,cl,B=100,balanced=FALSE,mat.samp=NULL,delta=(1:10)/5,med.fdr=TRUE,
		s0=NA,alpha.s0=seq(0,1,.05),include.s0=TRUE,p0=NA,lambda.p0=1,vec.lambda.p0=(0:95)/100,
		na.rm=FALSE,graphic.fdr=TRUE,thres.fdr=seq(0.5,2,0.5),ngenes=NA,iteration=3,
		initial.delta=c(.1,seq(.2,2,.2),4),rand=NA){
	if(any(delta<=0))
		stop("Delta must be larger than 0.")
	X.name<-match.call()$data
	xy.out<-xy.cal(cl)
	x<-xy.out$x
	y<-xy.out$y
	paired<-xy.out$paired
	tmp<-sam.old(data,x,y=y,paired=paired,mat.samp=mat.samp,B=B,balanced=balanced,
		na.rm=na.rm,s0=s0,alpha.s0=alpha.s0,include.s0=include.s0,p0=p0,lambda.p0=lambda.p0,
		vec.lambda.p0=vec.lambda.p0,delta.fdr=delta,med.fdr=med.fdr,graphic.fdr=graphic.fdr,
		thres.fdr=thres.fdr,ngenes=ngenes,iteration=iteration,initial.delta=initial.delta,
		rand=rand)
	invisible(structure(list(X.name=X.name,d=tmp$d,d.sort=tmp$d.sort,s=tmp$s,d.bar=tmp$d.bar,d.perm=tmp$d.perm,
		mat.samp=tmp$mat.samp,s0=tmp$s0,p0=tmp$p0,FDR=tmp$FDR,fdr.ngenes=tmp$fdr.ngenes,
		delta.ngenes=tmp$delta.ngenes,med.fdr=tmp$med.fdr,x=x,y=y,paired=paired)))
}