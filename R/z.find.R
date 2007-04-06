`z.find` <-
function(data,cl,B=100,var.equal=FALSE,B.more=0.1,B.max=30000){
	data<-as.matrix(data)
	mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("No missing values allowed.")
	adjust.out<-adjust.for.mt(data,cl,var.equal=var.equal,eb=TRUE)
	data<-adjust.out$X
	cl<-adjust.out$cl.mt
	msg<-adjust.out$msg
	type.mt<-adjust.out$type.mt
	z.fun<-switch(EXPR=type.mt,
		t=function(data,cl) computeRS(data,cl,"t"),
		t.equalvar=function(data,cl) computeRS(data,cl,"t.equalvar"),
		pairt=function(data,cl) computeRS(data,cl,"pairt"),
		f=function(data,cl) computeRS(data,cl,"f"))
	mt.out<-mt.teststat.num.denum(data,cl,test=type.mt)
	mat.samp<-setup.mat.samp(cl,type.mt,B=B,B.more=B.more,B.max=B.max)
	if(type.mt=="pairt"){
		le.cl<-length(cl)/2
		data<-data[,2*(1:le.cl)]-data[,2*(1:le.cl)-1]
		mat.samp<-mat.samp[,2*(1:le.cl)]
	}
	n.genes<-nrow(data)
	if(type.mt=="f"){
		n.classes<-length(unique(cl))
		z.norm<-qf(((1:n.genes)-0.5)/n.genes,n.classes-1,n.genes-n.classes)
	}
	else
		z.norm<-qnorm(((1:n.genes)-0.375)/(n.genes+0.25))
	return(list(r=mt.out$teststat.num,s=mt.out$teststat.denum,z.fun=z.fun,mat.samp=mat.samp,
		msg=msg,data=data,cl=cl,z.norm=z.norm,type.mt=type.mt))
}

