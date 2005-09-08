sam.dstat<-function(data,cl,var.equal=FALSE,B=100,med=FALSE,s0=NA,s.alpha=seq(0,1,.05),
		include.zero=TRUE,p0=NA,n.subset=10,mat.samp=NULL,B.more=0.1,B.max=30000,
		lambda=seq(0,0.95,0.05),ncs.value="max",ncs.weights=NULL,delta=NULL,n.delta=10,
		gene.names=dimnames(data)[[1]],q.version=1,R.fold=1,R.unlog=TRUE,na.replace=TRUE,
		na.method="mean",rand=NA){
	if(is(data,"exprSet")){
		require(affy)
		chip.name<-annotation(data)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
	}
	else
		chip.name<-""
	d.out<-d.stat(data,cl,var.equal=var.equal,B=B,s0=s0,s.alpha=s.alpha,include.zero=include.zero,
		mat.samp=mat.samp,B.more=B.more,B.max=B.max,med=med,n.subset=n.subset,
		gene.names=gene.names,R.fold=R.fold,R.unlog=R.unlog,na.replace=na.replace,
		na.method=na.method,rand=rand)
	if(is.na(p0))
		p0<-pi0.est(na.exclude(d.out$p),lambda=lambda,ncs.value=ncs.value,
			ncs.weights=ncs.weights)$p0
	mat.fdr<-stats.cal(d.out$d,d.out$d.bar,d.out$vec.false,p0,delta=delta,le.delta=n.delta)
	if(q.version%in%c(1,2))
		q.value<-qvalue.cal(d.out$p.value,p0,version=q.version)
	else
		q.value<-numeric(0)
	new("SAM",d=d.out$d,d.bar=d.out$d.bar,vec.false=d.out$vec.false,p.value=d.out$p.value,
		s=d.out$s,s0=d.out$s0,mat.samp=d.out$mat.samp,p0=p0,mat.fdr=mat.fdr,q.value=q.value,
		fold=d.out$fold,msg=d.out$msg,chip=chip.name)
}	 
 

