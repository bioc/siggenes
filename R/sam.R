sam <- function(data, cl, method=d.stat, control=samControl(), 
		gene.names=dimnames(data)[[1]], ...){
	if(is(data,"ExpressionSet")){
		require(affy)
		chip.name<-annotation(data)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
	}
	else
		chip.name<-""
	FUN<-match.fun(method)
	if(missing(cl)){
		usedFun <- deparse(substitute(method))
		if(usedFun %in% c("chisq.stat", "trend.stat")){
			if(!inherits(data, "list"))
				stop("data must be a list.")
			d.out <- FUN(data, ...)
		}
		else
			stop("cl needs to be specified.")
	}
	else{				
		if(is.factor(cl))
			cl<-as.character(cl)
		d.out<-FUN(data,cl,...)
	}
	if(is.na(control$p0))
		p0 <- pi0.est(na.exclude(d.out$p), lambda=control$lambda,
			ncs.value=control$ncs.value, ncs.weights=control$ncs.weights)$p0
	else
		p0 <- control$p0
	mat.fdr <- stats.cal(d.out$d, d.out$d.bar, d.out$vec.false, p0, 
		delta=control$delta, le.delta=control$n.delta)
	if(!is.null(gene.names)){
		names(d.out$d)<-names(d.out$p.value)<-names(d.out$vec.false)<-substring(gene.names,1,50)
		if(length(d.out$s)==length(d.out$d))
			names(d.out$s)<-names(d.out$d)
	}
	if(control$q.version %in% c(1,2))
		q.value <- qvalue.cal(d.out$p.value, p0, version=control$q.version)
	else
		q.value <- numeric(0)
	new("SAM",d=d.out$d,d.bar=d.out$d.bar,vec.false=d.out$vec.false,p.value=d.out$p.value,
		s=d.out$s,s0=d.out$s0,mat.samp=d.out$mat.samp,p0=p0,mat.fdr=mat.fdr,q.value=q.value,
		fold=d.out$fold,msg=d.out$msg,chip=chip.name)
} 
