fuzzy.stat <- function(data, cl, type=c("asymptotic", "permutation", "abf"), W=NULL, logbase=exp(1),
		addOne=TRUE, B=100, B.more=0.1, B.max=30000, n.subset=10, rand=NA){
	require(scrime, quietly=TRUE)
	if(is.data.frame(data))
		data <- as.matrix(data)
	if(!is.matrix(data))
		stop("data must be either a matrix or a data frame.")
	if(!is.na(logbase) && logbase <= 1)
		stop("logbase must be larger than 1.")
	out <- rowTrendFuzzy(y=cl, mat.fuzzy=data)
	type <- match.arg(type)
	if(type=="abf"){
		if(is.null(W)){
			warning("Since W is missing, it is set to (log(2)/qnorm(0.95))^2.")
			W <- (log(2) / qnorm(0.95))^2
		}
		z <- abf(out$theta, out$varTheta, W=W, numerator=1)
		if(!is.na(logbase))
			z <- log(z+addOne, base=logbase)
		sd1 <- sqrt(out$varTheta + W)
	}
	else{
		z <- out$stat
		sd1 <- NULL
	}
	s <- sqrt(out$varTheta)
	if(type=="asymptotic"){
		null.out <- fuzzy.null.approx(out$rawp)
		mat.samp <- matrix(numeric(0))
	}
	else{
		mat.samp <- setup.mat.samp(cl, "t", B=B, B.more=B.more, B.max=B.max, rand=rand)
		null.out <- fuzzy.null(data, mat.samp, z, s, sd1, n.subset=n.subset, logbase=logbase, add1=addOne)	 
	}
	msg <- paste("SAM for Fuzzy Genotype Calls", if(type=="abf") "(Based on Approximate Bayes Factors)", "\n\n")
	list(d=z, d.bar=null.out$d.bar, p.value=null.out$p.value, vec.false=null.out$vec.false, 
		s=s, s0=numeric(0), mat.samp=mat.samp, msg=msg, fold=numeric(0))
}

fuzzy.null.approx <- function(rawp){
	m <- length(na.exclude(rawp))
	d.bar <- qnorm(((1:m) - 0.5) / m)
	vec.false <- m * rawp / 2
	list(d.bar=d.bar, p.value=rawp, vec.false=vec.false)
}

fuzzy.null <- function(mat.fuzzy, mat.samp, d, s, sd1, n.subset=10, logbase=exp(1), add1=TRUE){
	if(!is.null(sd1))
		return(fuzzy.null.abf(mat.fuzzy, mat.samp, d, s, sd1, n.subset=n.subset, logbase=logbase, add1=add1))
	B <- nrow(mat.samp)
	vecB <- unique(c(seq(1, B, n.subset), B+1))
	n.int <- length(vecB) - 1
	n.var <- nrow(mat.fuzzy)
	mat.samp <- mat.samp - mean(mat.samp[1,])
	dpos <- d >= 0
	d <- abs(d)
	d.rank <- rank(-d, ties.method="first")
	vec.dperm <- numeric(n.var)
	vec.pos <- vec.neg <- -n.int * d.rank
	for(i in 1:n.int){
		tmpmat <- t(mat.samp[vecB[i]:(vecB[i+1] - 1),, drop=FALSE])
		mat.dperm <- mat.fuzzy %*% tmpmat * s
		mat.dperm <- apply(mat.dperm, 2, sort)
		vec.dperm <- vec.dperm + rowSums(mat.dperm, na.rm=TRUE)
		pos <- mat.dperm >= 0
		tmp <- rank(-c(mat.dperm[pos], d), na.last=NA, ties.method="first")
		vec.pos <- vec.pos + tmp[length(tmp) - ((n.var-1):0)] 
		tmp <- rank(c(mat.dperm[!pos], -d), na.last=NA, ties.method="first")
		vec.neg <- vec.neg + tmp[length(tmp) - ((n.var-1):0)] 
	}
	d.bar <- vec.dperm / B
	p.value <- (vec.pos + vec.neg) / (n.var * B)
	vec.false <- numeric(n.var)
	vec.false[dpos] <- vec.pos[dpos]
	vec.false[!dpos] <- vec.neg[!dpos]
	structure(list(d.bar=d.bar, vec.false=vec.false/B, p.value=p.value))
}
		
fuzzy.null.abf <- function(mat.fuzzy, mat.samp, d, s, sd1, n.subset=10, logbase=exp(1), add1=TRUE){
	B <- nrow(mat.samp)
	vecB <- unique(c(seq(1, B, n.subset), B+1))
	n.int <- length(vecB) - 1
	n.var <- nrow(mat.fuzzy)
	mat.samp <- mat.samp - mean(mat.samp[1,])
	vec.dperm <- numeric(n.var)
	d.rank <- rank(-d, ties.method="first")
	vec.false <- -n.int * d.rank
	var0 <- s * s
	for(i in 1:n.int){
		tmpmat <- mat.fuzzy %*% t(mat.samp[vecB[i]:(vecB[i+1] - 1),, drop=FALSE]) * var0 
		mat.dperm <- dnorm(tmpmat, sd=sd1) / dnorm(tmpmat, sd=s)
		if(!is.na(logbase))
			mat.dperm <- log(mat.dperm + add1, base=logbase)
		mat.dperm <- apply(mat.dperm, 2, sort)
		vec.dperm <- vec.dperm + rowSums(mat.dperm)
		tmp <- rank(-c(as.vector(mat.dperm), d), ties.method="first")
		vec.false <- vec.false + tmp[length(tmp)-((n.var-1):0)]
	}
	d.bar <- vec.dperm / B
	vec.false <- vec.false / B
	p.value <- vec.false / n.var
	structure(list(d.bar=d.bar, vec.false=vec.false, p.value=p.value))
}


		





