trend.stat <- function(data, cl, ...) UseMethod("trend.stat")

trend.stat.data.frame <- function(data, cl, ...){
	data <- as.matrix(data)
	trend.stat(data, cl, ...)
}

trend.stat.default <- function(data, cl, catt=TRUE, approx=TRUE, B=100, B.more=0.1, B.max=50000,
		n.subset=10, rand=NA, ...){
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
	if(approx){
		null.out <- cat.null.approx(stats, 2, 2)
		mat.samp <- matrix(numeric(0))
	}
	else{
		mat.samp <- setupMatSamp2(cl, n.cat-1, B=B, B.more=B.more, B.max=B.max, 
			rand=rand)
		facObs <- if(catt) (n.obs-1)^2 / n.obs else n.obs 
		null.out <- trend.null(data, mat.samp, stats, facObs, n.subset, B) 
	}
	msg <- c("SAM Analysis of Linear Trend\n\n",
		paste("Test Statistic: ", 
		ifelse(catt, "Cochran-Armitage", "MSquare"), "\n", sep=""),
		paste("Null Distribution: ",
		ifelse(approx, "Approximation by ChiSquare Distribution with 1 df", 
		paste("Estimated Based on", B, "Permutations")), "\n\n", sep=""))
	structure(list(d=stats, d.bar=null.out$d.bar, p.value=null.out$p.value,
		vec.false=null.out$vec.false, discrete=FALSE, s=numeric(0),
		s0=numeric(0), mat.samp=mat.samp, msg=msg, fold=numeric(0)))
}
		

setupMatSamp2 <- function(cl, nFrom0, B=100, B.more=0.1, B.max=50000, rand=NA){
	uni.cl <- unique.default(cl)
	tab <- rep(0:nFrom0, tabulate(factor(cl, levels=uni.cl)))
	tmp <- setup.mat.samp(tab, "f", B=B, B.more=B.more, B.max=B.max, rand=rand) + 1
	matrix(uni.cl[tmp], B)
} 


trend.null <- function(x, mat.perm, stats, facObs, n.subset, B){
	vecB <- c(seq(1, B, n.subset), B + 1)
	vecB <- unique.default(vecB)
	n.B <- length(vecB) - 1
	n.var <- length(stats)
	vec.dperm <- vec.false <- numeric(n.var)
	d.rank <- rank(-stats, ties="first")
	for(i in 1:n.B){
		tmp <- compPermTrendStat(x, mat.perm[vecB[i]:(vecB[i+1]-1),, drop=FALSE], facObs)
		vec.dperm <- vec.dperm + rowSums(tmp)
		tmp2 <- c(as.vector(tmp), stats)
		vec.false <- vec.false + rank(-tmp2, ties="first")[length(tmp) + (1:n.var)] - d.rank
	}
	d.bar <- vec.dperm / B
	vec.false <- vec.false / B
	p.value <- vec.false / n.var
	return(list(d.bar=d.bar, vec.false=vec.false, p.value=p.value))
}


compPermTrendStat <- function(x, perms, facObs){
	tmp <- x %*% t(perms)
	tmp <- tmp * tmp / facObs
	apply(tmp, 2, sort)
}
	


trend.stat.list <- function(data, cl, catt=TRUE, approx=TRUE, B=100, B.more=0.1, B.max=50000, 
		n.subset=10,rand=NA, ...){
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
	out <- if(catt) rowCATTs(data[[1]], data[[2]]) 
		else rowMsquares(listTables=data, clScores=cl)
	n.var <- length(out$stats)
	d.bar <- qchisq(((1:n.var) - 0.5) / n.var, 1)
	vec.false <- out$rawp * n.var
	msg <- c("SAM Analysis of Linear Trend\n\n",
		paste("Test Statistic:", ifelse(catt, "Cochran-Armitage", "MSquare"), "\n"),
		"Null Distribution: Approximation by ChiSquare Distribution with 1 df\n\n")
	structure(list(d=out$stats, d.bar=d.bar, p.value=out$rawp, vec.false=vec.false,
		discrete=FALSE, s=numeric(0), s0=numeric(0), mat.samp=matrix(numeric(0)),
		msg=msg, fold=numeric(0)))
}


		


