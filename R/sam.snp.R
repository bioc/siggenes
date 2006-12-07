sam.snp<-function(data,cl,B=100,approx=FALSE,delta=NULL,n.delta=10,p0=NA,lambda=seq(0,0.95,0.05),
		ncs.value="max",ncs.weights=NULL,gene.names=dimnames(data)[[1]],q.version=1,
		check.levels=TRUE,check.for.NN=FALSE,lev=NULL,B.more=0.1,B.max=50000,
		n.subset=10,rand=NA){
	d.out<-cat.stat(data,cl,B=B,check.levels=check.levels,check.for.NN=check.for.NN,
		approx=approx,lev=lev,B.more=B.more,B.max=B.max,n.subset=n.subset,rand=rand)
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
		fold=d.out$fold,msg=d.out$msg,chip="")
} 
 

