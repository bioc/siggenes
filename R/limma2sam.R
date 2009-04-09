limma2sam <- function(fit, coef, moderate = TRUE, sam.control=samControl()){
	if(!is(fit,"MArrayLM"))
		stop("fit must be an object of class MArrayLM.")
	if(length(coef)!=1)
		stop("coef must be of length 1.")
	if(moderate){
		if(is.null(fit$t))
			stop("fit must be the output of eBayes if moderate = TRUE.")
		d <- fit$t[,coef]
		df <- fit$df.prior + fit$df.residual
		p.value <- fit$p.value[,coef]
	}
	else{
		d <- fit$coefficients[,coef] / fit$stdev.unscaled[,coef] / fit$sigma
		df <- fit$df.residual
		p.value <- 2 * pt(-abs(d), df=df)
	}
	m <- sum(!is.na(d))
	d.bar <- qt(ppoints(m), df=df)
	vec.false <- m * p.value / 2
	num <- fit$coefficients[,coef]
	s <- num / d
	p0 <- sam.control$p0
	q.version <- sam.control$q.version
	if(is.na(p0))
		p0 <- pi0.est(na.exclude(p.value), lambda=sam.control$lambda, 
			ncs.value=sam.control$ncs.value, 
			ncs.weights=sam.control$ncs.weigths)$p0
	if(p0<=0 | p0>1)
		stop("p0 must be between 0 and 1.")
	mat.fdr <- stats.cal(d, d.bar, vec.false, p0, delta=sam.control$delta, 
		le.delta=sam.control$n.delta)
	if(q.version %in% (1:2))
		q.value <- qvalue.cal(p.value, p0, version=q.version)
	else
		q.value <- numeric(0)
	msg <- c("SAM Analysis for LIMMA Object\n\n", paste("Column: ", coef, "\n", 
		"T-Statistic: ", ifelse(moderate, "Moderated", "Ordinary"), "\n\n", sep=""))
	new("SAM", d=d, d.bar=d.bar, vec.false=vec.false, p.value=p.value, s=s, s0=numeric(0),
		mat.samp=matrix(numeric(0)), p0=p0, mat.fdr=mat.fdr, q.value=q.value,
		fold=2^num, msg=msg, chip="")
}

samControl <- function(delta=NULL, n.delta=10, p0=NA, lambda=seq(0,0.95,0.05), ncs.value="max",
		ncs.weights=NULL, q.version=1){
	list(delta=delta, n.delta=n.delta, p0=p0, lambda=lambda, ncs.value=ncs.value,
		ncs.weights=ncs.weights, q.version=q.version)
}

ebamControl <- function(p0=NA, p0.estimation=c("splines","interval","adhoc"), lambda=NULL,
		ncs.value="max", use.weights=FALSE){
	list(p0=p0, p0.estimation=match.arg(p0.estimation), lambda=lambda, ncs.value=ncs.value,
		use.weights=use.weights)
}

limma2ebam <- function(fit, coef, moderate=TRUE, delta=0.9, ebam.control=ebamControl()){
	if(!is(fit,"MArrayLM"))
		stop("fit must be an object of class MArrayLM.")
	if(length(coef)!=1)
		stop("coef must be of length 1.")
	if(moderate){
		if(is.null(fit$t))
			stop("fit must be the output of eBayes if moderate = TRUE.")
		z <- fit$t[,coef]
		df <- fit$df.prior + fit$df.residual
		p.value <- fit$p.value[,coef]
	}
	else{
		z <- fit$coefficients[,coef] / fit$stdev.unscaled[,coef] / fit$sigma
		df <- fit$df.residual
		p.value <- 2 * pt(-abs(z), df=df)
	}
	if(any(is.na(z)))
		stop("No missing z-values allowed.")
	m <- length(z)
	z.dens <- denspr(z)$y
	vec.pos <- m * p.value / 2
	z.null <- dt(z, df=df)
	ratio <- z.null / z.dens
	msg <- paste("EBAM Analysis for LIMMA Object\n\n", "Column: ", coef, "\n", 
		"T-Statistic: ", ifelse(moderate, "Moderated", "Ordinary"), "\n\n", sep="")
	p0 <- ebam.control$p0
	type.p0 <- ebam.control$p0.estimation
	lambda <- ebam.control$lambda
	if(is.na(p0)){
		if(is.null(lambda))
			lambda <- if(type.p0=="splines") seq(0, 0.95, 0.05)  else 0.5
		w <- if(ebam.control$use.weights) 1 - lambda  else NULL
		p0 <- if(type.p0=="adhoc") min(1/ratio)
			else pi0.est(p.value, lambda=lambda, ncs.value=ebam.control$ncs.value,
				ncs.weights=w)$p0
	}
	if(p0<=0 | p0>1)
		stop("p0 must be between 0 and 1.")
	posterior <- 1 - p0 * ratio
	if(any(posterior<0)){
		warning("Some of the posterior probabilities are less than zero.\n",
			"These probabilities are set to zero.")
		posterior[posterior<0] <- 0
	}
	mat <- compNumber(z, posterior, p0, NA, delta=delta, vec.pos=vec.pos, vec.neg=vec.pos)
	new("EBAM", z=z, posterior=posterior, p0=p0, local=1-posterior, mat.fdr=mat,
		a0=numeric(0), mat.samp=matrix(numeric(0)), vec.pos=vec.pos, vec.neg=vec.pos,
		msg=msg, chip="")
}

