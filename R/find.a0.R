`find.a0` <- function(data, cl, method=z.find, B=100, delta=0.9, quan.a0=(0:5)/5,
		include.zero=TRUE, control=find.a0Control(), gene.names=dimnames(data)[[1]],
		rand=NA,...){
	if(length(delta)>1)
		stop("For find.a0, only one value of delta is allowed.\n",
			"Use print to obtain the number of interesting genes for other values ",
			"of delta.")
	if(is(data,"ExpressionSet")){
		# requireNamespace("affy",quietly=TRUE)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		chip.name <- data@annotation
		data<-exprs(data)
	}
	else
		chip.name<-""
	if(length(cl)!=ncol(data))
		stop("The number of columns of data must be equal to the length of cl.")
	FUN<-match.fun(method)
	if(!is.na(rand))
		set.seed(rand)
	out<-FUN(data,cl,B=B,...)
	tmp.names<-names(out)
	if(any(!c("r","s","z.fun","mat.samp")%in%tmp.names))
		stop("The function specified by method must return a list consisting of\n",
			"objects called r, s, z.fun, and mat.samp.")
	r<-out$r
	s<-out$s
	z.fun<-out$z.fun
	mat.samp<-out$mat.samp
	n.genes<-length(r)
	z.norm<-if(is.null(out$z.norm)) qnorm(((1:n.genes)-0.375)/(n.genes+0.25))
		else out$z.norm
	if(substitute(method)=="z.find"){
		data<-out$data
		cl<-out$cl
	}
	msg<-if(any(names(out)=="msg")) out$msg  else "" 
	rm(out)
	vec.a0<-checkQuantiles(s,quan.a0=quan.a0,include.zero=include.zero)
	z.mat<-r/outer(s,vec.a0,"+")
	succ.out <- apply(z.mat, 2, getSuccesses, n.interval=control$n.interval)
	tmp<-lapply(succ.out,function(x) x$success)
	mat.success <- matrix(unlist(tmp), nrow=control$n.interval)
	ints<-lapply(succ.out,function(x) x$interval)
	tmp<-lapply(succ.out,function(x) x$center)
	mat.center <- matrix(unlist(tmp), nrow=control$n.interval)
	succ.norm <- getSuccesses(z.norm, n.interval=control$n.interval)
	fail.mats <- compFailureMat2(data, mat.samp, z.mat, ints, z.norm, 
		succ.norm$interval, z.fun, vec.a0, n.chunk=control$n.chunk,
		n.interval=control$n.interval)
	mat.fn<-fail.mats$mat.fail.norm
	n.a0<-length(vec.a0)
	mat.ratio<-matrix(0,n.genes,n.a0)
	df.ratio <- control$df.ratio
	if(is.null(df.ratio))
		df.ratio<-ifelse(any(r<0),5,3)
	for(i in 1:n.a0)
		mat.ratio[,i]<-compRatio(succ.norm$center,succ.norm$success,mat.fn[,i],
			df=df.ratio,z=z.norm)$ratio
	type.p0 <- control$p0.estimation
	if(type.p0=="adhoc")
		vec.p0<-apply(mat.ratio,2,function(x) min(1/x))
	else
		vec.p0 <- apply(mat.fn, 2, pi0.est2, success=succ.norm$success,
			z=succ.norm$center, type=type.p0, lambda=control$lambda,
			ncs.value=control$ncs.value, use.weights=control$use.weights)
	mat.posterior<-1-rep(vec.p0,e=n.genes)*mat.ratio
	mat.posterior[mat.posterior<0]<-0
	names(vec.a0)<-c(if(vec.a0[1]==0) "-", quan.a0)
	a0.out<-makeA0mat(z.norm,mat.posterior,vec.p0,vec.a0,B,delta=delta)
	colnames(z.mat)<-colnames(mat.posterior)<-colnames(mat.success)<-
		colnames(mat.center)<-names(vec.p0)<-round(vec.a0,4)
	if(!is.null(gene.names))
		rownames(z.mat)<-rownames(mat.posterior)<-gene.names
	new("FindA0",mat.z=z.mat,mat.posterior=mat.posterior,mat.center=mat.center,
		mat.success=mat.success,mat.failure=fail.mats$mat.failure,
		z.norm=z.norm,p0=vec.p0,mat.a0=a0.out$tab,
		mat.samp=mat.samp,vec.a0=vec.a0,suggested=a0.out$suggest,delta=delta,
		df.ratio=df.ratio,msg=msg,chip=chip.name)
}

