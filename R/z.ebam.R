`z.ebam` <-
function(data,cl,a0=NULL,quan.a0=NULL,B=100,delta=.9,var.equal=FALSE,B.more=0.1,B.max=30000,
		gene.names=dimnames(data)[[1]],n.subset=10,fast=FALSE,n.interval=139,
		df.ratio=NULL,rand=NA){
	if(!is.na(rand))
		set.seed(rand)
	out<-z.find(data,cl,B=B,var.equal=var.equal,B.more=B.more,B.max=B.max)
	a0<-checkA0(out$s,a0=a0,quan.a0=quan.a0)
	z<-out$r/(a0+out$s)
	data<-out$data
	mat.samp<-out$mat.samp
	msg<-out$msg
	type.mt<-out$type.mt
	out<-getSuccesses(z,n.interval=n.interval)
	fail.out<-compFailure(data,mat.samp,z,out$interval,a0,type.mt,n.subset=n.subset,fast=fast)
	if(is.null(df.ratio))
		df.ratio<-ifelse(any(z<0),5,3)
	ratio<-compRatio(out$center,out$success,fail.out$vec.fail,df=df.ratio,z=z)$ratio
	if(fast)
		return(list(z=z,ratio=ratio,a0=a0,success=out$success,failure=fail.out$vec.fail,
			center=out$center,mat.samp=mat.samp,msg=msg))
	else
		return(list(z=z,ratio=ratio,a0=a0,vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B,mat.samp=mat.samp,msg=msg))		
}

