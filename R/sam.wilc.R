sam.wilc<-function(data,cl,delta=NULL,n.delta=10,p0=NA,lambda=seq(0,0.95,0.05),ncs.value="max",
		ncs.weights=NULL,gene.names=dimnames(data)[[1]],q.version=1,R.fold=1,R.unlog=TRUE,
		na.replace=TRUE,na.method="mean",approx50=TRUE,check.ties=FALSE,rand=NA){
	if(is(data,"exprSet")){
		require(affy)
		chip.name<-annotation(data)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
	}
	else
		chip.name<-""
	d.out<-wilc.stat(data,cl,gene.names=gene.names,R.fold=R.fold,R.unlog=R.unlog,
		na.replace=na.replace,na.method=na.method,approx50=approx50,check.ties=check.ties,
		rand=rand)
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
		 
 

