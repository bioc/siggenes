trend.ebam <- function(data, cl, ...) UseMethod("trend.ebam")

trend.ebam.data.frame <- function(data, cl, ...){
	data <- as.matrix(data)
	trend.ebam(data, cl, ...)
}

trend.ebam.default <- function(data, cl, catt=TRUE, approx=TRUE, n.interval=NULL, df.dens=NULL,
		knots.mode=NULL, type.nclass="wand",B=100, B.more=0.1, B.max=50000,
		n.subset=10, fast=FALSE, df.ratio=3, rand=NA, ...){
	if(!is.matrix(data) || !is.numeric(data))
		stop("data must be a numeric matrix or data frame.")
	if(missing(cl))
		stop("cl must be specified if data is a matrix or data frame.")
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.")
	if(!is.numeric(data))
		stop("data must be numeric.")
	if(!is.numeric(cl))
		stop("cl must be numeric.")
	if(length(cl) != ncol(data))
		stop("The number of columns of data must be equal to the length of cl.")
	n.cat <- length(unique(cl))
	if(n.cat < 2)
		stop("cl must consist of at least two different scores/levels.")
	if(n.cat > 10)
		stop("cl consists of more than ten different scores.")
	if(any(is.na(data))){
		if(!approx)
			stop("No missing values allowed when approx = FALSE.")
		n.obs <- rowSums(!is.na(data)) - 1
	}
	else
		n.obs <- ncol(data) - 1
	if(any(n.obs <= 1))
		stop("Each row of data must contain at least two non-missing values.")
	cl <- (cl - mean(cl)) / sd(cl) 
	data <- data - rowMeans(data, na.rm=TRUE)
	rSd <- rowSums(data * data, na.rm = TRUE) / n.obs
	data <- data / sqrt(rSd)
	stats <- as.vector(data %*% cl)
	stats <- stats / n.obs
	catt <- catt & (n.cat==2)
	if(catt)
		n.obs <- n.obs + 1
	stats <- stats * stats * n.obs
	msg <- c("EBAM Analysis of Linear Trend\n\n",
		paste("Test Statistic: ", ifelse(catt, "Cochran-Armitage", "MSquare"), "\n",
			sep=""),
		paste("Null Distribution: ",
			ifelse(approx, "Approximation by ChiSquare Distribution with 1 df",
			paste("Estimated Based on", B, "Permutations")), "\n\n", sep=""))
	if(approx){
		fail.out <- cat.null.approx2(stats, 2, 2, n.interval=n.interval, 
			df.dens=df.dens, knots.mode=knots.mode, type.nclass=type.nclass)
		return(list(z=stats, ratio=fail.out$ratio, vec.pos=fail.out$vec.pos,
			vec.neg=fail.out$vec.neg, msg=msg))
	}
	if(is.null(n.interval))
		n.interval <- 139
	out <- getSuccesses(stats, n.interval=n.interval)
	mat.samp <- setupMatSamp2(cl, n.cat-1, B=B, B.more=B.more, B.max=B.max, rand=rand)
	B <- nrow(mat.samp)
	facObs <- if(catt) (n.obs)^2 / n.obs   else n.obs
	fail.out <- compFailure(data, mat.samp, stats, out$interval, n.subset=n.subset,
		fast=fast, n.cat=-facObs)
	ratio <- compRatio(out$center, out$success, fail.out$vec.fail, df=df.ratio,
		z=stats)$ratio
	if(fast)
		return(list(z=stats, ratio=ratio, success=out$success, failure=fail.out$vec.fail,
			center=out$center, mat.samp=mat.samp,msg=msg))
	else
		return(list(z=stats, ratio=ratio, vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B, mat.samp=mat.samp, msg=msg))
}

trend.ebam.list <- function(data, cl, catt=TRUE, approx=TRUE, n.interval=NULL, df.dens=NULL, 
		knots.mode=NULL, type.nclass="wand", ...){
	require(scrime)
	if(missing(cl))
		cl <- NULL
	n.cat <- length(data)
	if(!is.null(cl) && length(cl) != n.cat)
		stop("If data is a list, then cl (if specified) must contain one score\n",
			"for each entry of data.")
	if(!approx)
		stop("Currently only approx=TRUE available if data is a list.")
	catt <- catt & n.cat==2
	stats <- if(catt) rowCATTs(data[[1]], data[[2]], add.pval=FALSE)
		else rowMsquares(listTables=data, clScores=cl, add.pval=FALSE)
	out <- cat.null.approx2(stats, 2, 2, n.interval=n.interval, df.dens=df.dens,
		knots.mode=knots.mode, type.nclass=type.nclass)
	msg <- c("EBAM Analysis of Linear Trend\n\n",
		paste("Test Statistic: ", ifelse(catt, "Cochran-Armitage", "MSquare"), "\n",
			sep=""),
		"Null Distribution: Approximation by ChiSquare Distribution with 1 df\n\n")
	structure(list(z=stats, ratio=out$ratio, vec.pos=out$vec.pos, vec.neg=out$vec.neg,
		msg=msg))
}
	