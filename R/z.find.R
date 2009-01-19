`z.find` <-
function(x,y,B=100,var.equal=FALSE,B.more=0.1,B.max=30000){
	x<-as.matrix(x)
	mode(x)<-"numeric"
	if(any(is.na(x)))
		stop("No missing values allowed.")
	adjust.out<-adjust.for.mt(x,y,var.equal=var.equal,eb=TRUE)
	x<-adjust.out$X
	y<-adjust.out$cl.mt
	msg<-adjust.out$msg
	type.mt<-adjust.out$type.mt
	z.fun<-switch(EXPR=type.mt,
		t=function(x,y) computeRS(x,y,"t"),
		t.equalvar=function(x,y) computeRS(x,y,"t.equalvar"),
		pairt=function(x,y) computeRS(x,y,"pairt"),
		f=function(x,y) computeRS(x,y,"f"))
	mt.out<-mt.teststat.num.denum(x,y,test=type.mt)
	mat.samp<-setup.mat.samp(y,type.mt,B=B,B.more=B.more,B.max=B.max)
	if(type.mt=="pairt"){
		le.cl<-length(y)/2
		x<-x[,2*(1:le.cl)]-x[,2*(1:le.cl)-1]
		mat.samp<-mat.samp[,2*(1:le.cl)]
	}
	n.genes<-nrow(x)
	if(type.mt=="f"){
		n.classes<-length(unique(y))
		z.norm<-qf(((1:n.genes)-0.5)/n.genes,n.classes-1,n.genes-n.classes)
	}
	else
		z.norm<-qnorm(((1:n.genes)-0.375)/(n.genes+0.25))
	return(list(r=mt.out$teststat.num,s=mt.out$teststat.denum,z.fun=z.fun,mat.samp=mat.samp,
		msg=msg,x=x,y=y,z.norm=z.norm,type.mt=type.mt))
}

