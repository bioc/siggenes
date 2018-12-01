fuzzy.ebam <- function(data, cl, type=c("asymptotic", "permutation", "abf"), W=NULL,
		logbase=exp(1), addOne=TRUE, df.ratio=NULL, n.interval=NULL, df.dens=5, 
		knots.mode=TRUE, type.nclass=c("FD","wand", "scott"), fast=FALSE,
		B=100, B.more=0.1, B.max=30000, n.subset=10, rand=NA){
	# requireNamespace("scrime", quietly=TRUE)
	if(is.data.frame(data))
		data <- as.matrix(data)
	if(!is.matrix(data))
		stop("data must be either a matrix or a data frame.")
	if(!is.na(logbase) && logbase <= 1)
		stop("logbase must be larger than 1.")
	if(!is.logical(addOne))
		stop("addOne must be a logic value.")
	type <- match.arg(type)
	type.nclass <- match.arg(type.nclass)
	out <- rowTrendFuzzy(y=cl, mat.fuzzy=data)
	msg <- paste("EBAM for Fuzzy Genotype Calls", 
		if(type=="abf") "(Based on Approximate Bayes Factors)", "\n\n")
	if(type=="asymptotic"){
		z.dens <- denspr(out$stat, n.interval=n.interval, df=df.dens,
			knots.mode=knots.mode, type.nclass=type.nclass)$y
		vec.pos <- nrow(data) * out$rawp/2
		z.null <- dnorm(out$stat)
		return(list(z=out$stat, ratio=z.null/z.dens, vec.pos=vec.pos, vec.neg=vec.pos,
			msg=msg))
	}
	if(type=="abf"){
		if(is.null(W)){
			warning("Since W is missing, it is set to (log(2)/qnorm(0.95))^2.")
			W <- (log(2) / qnorm(0.95))^2
		}
		z <- abf(out$theta, out$varTheta, W=W, numerator=1)
		if(!is.na(logbase))
			z <- log(z + addOne, base=logbase)
		sd1 <- sqrt(out$varTheta + W)
	}
	else{
		z <- out$stat
		sd1 <- NULL
	}
	if(is.null(df.ratio))
		df.ratio <- ifelse(type=="abf", 3, 5)
	if(is.null(n.interval)){
		tmpfun <- match.fun(paste("nclass", type.nclass, sep="."))
		n.interval <- max(tmpfun(z), 139)
	}
	if(n.interval < 20)
		stop("n.interval should be at least 20.")
	succ.out <- getSuccesses(z, n.interval=n.interval)
	mat.samp <- setup.mat.samp(cl, "t", B=B, B.more=B.more, B.max=B.max, rand=rand)
	B <- nrow(mat.samp)
	fail.out <- compFailFuzzy(data, mat.samp, z, sqrt(out$varTheta), succ.out$interval,
		n.subset=n.subset, fast=fast, sd1=sd1, logbase=logbase, add1=addOne)
	ratio <- compRatio(succ.out$center, succ.out$success, fail.out$vec.fail,
		df=df.ratio, z=z)$ratio
	if(fast)
		return(list(z=z, ratio=ratio, success=succ.out$success, failure=fail.out$vec.fail,
			center=succ.out$center, mat.samp=mat.samp, msg=msg))
	structure(list(z=z, ratio=ratio, vec.pos=fail.out$vec.pos/B, vec.neg=fail.out$vec.neg/B,
		mat.samp=mat.samp, msg=msg))
}


compFailFuzzy <- function(X, mat.samp, z, s, interval, n.subset=10, fast=FALSE, sd1=NULL, logbase=exp(1),
		add1=TRUE){
	if(!is.null(sd1))
		return(compFailABF(X, mat.samp, z, s, sd1, interval, n.subset=n.subset, fast=fast, 
			logbase=logbase, add1=add1))
	B <- nrow(mat.samp)
	seq.samp <- unique(c(seq(1, B, ceiling(B/n.subset)), B+1))
	n.int <- length(interval) - 1
	n.seq <- length(seq.samp) - 1
	vec.fail <- numeric(n.int)
	n.row <- nrow(X)
	z.range <- range(z)
	mat.samp <- mat.samp - mean(mat.samp[1,])
	if(!fast){
		vec.pos <- vec.neg <- numeric(n.row)
		z.rank <- rank(-abs(z), ties.method="first")
	}
	else
		vec.pos <- vec.neg <- NULL
	for(i in 1:n.seq){
		tmpmat <- t(mat.samp[seq.samp[i]:(seq.samp[i+1] - 1),, drop=FALSE])
		z.perm <- as.vector(X %*% tmpmat * s)
		tmp <- getFailure(z.perm, z, interval, z.range=z.range, n.interval=n.int)
		vec.fail <- vec.fail + tmp
		if(!fast){
			tmp <- compFalse(z, z.perm, z.rank, n.row)
			vec.pos <- vec.pos + tmp$vec.pos
			vec.neg <- vec.neg + tmp$vec.neg
		}
	}
	structure(list(vec.fail=vec.fail, vec.pos=vec.pos, vec.neg=vec.neg))
}


compFailABF <- function(X, mat.samp, z, s, sd1, interval, n.subset=10, fast=FALSE, logbase=exp(1), add1=TRUE){
	B <- nrow(mat.samp)
	seq.samp <- unique(c(seq(1, B, ceiling(B/n.subset)), B+1))
	n.int <- length(interval) - 1
	n.seq <- length(seq.samp) - 1
	n.row <- nrow(X)
	mat.samp <- mat.samp - mean(mat.samp[1,])
	z.range <- range(z)
	var0 <- s*s
	vec.fail <- numeric(n.int)
	if(!fast){
		vec.neg <- numeric(n.row)
		vec.pos <- -n.seq * rank(-z, ties.method="first")
	}
	else
		vec.pos <- vec.neg <- NULL
	for(i in 1:n.seq){
		tmpvec <- as.vector(X %*% t(mat.samp[seq.samp[i]:(seq.samp[i+1] - 1),, drop=FALSE]) * var0)
		z.perm <- dnorm(tmpvec, sd=sd1) / dnorm(tmpvec, sd=s)
		if(!is.na(logbase))
			z.perm <- log(z.perm + add1, base=logbase)
		tmp <- getFailure(z.perm, z, interval, z.range=z.range, n.interval=n.int)
		vec.fail <- vec.fail + tmp
		if(!fast)
			vec.pos <- vec.pos + rank(-c(z.perm, z), ties.method="first")[length(z.perm) + (1:n.row)]
	}
	structure(list(vec.fail=vec.fail, vec.pos=vec.pos, vec.neg=vec.neg))
}




		