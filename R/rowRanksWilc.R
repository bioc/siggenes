rowRanksWilc <- function(X, cl){
	if(length(cl)!=ncol(X))
		stop("The length of cl must be equal to the number of columns of X.")
	ids <- which(cl == 1)
	Xr <- matrix(0, nrow(X), length(ids))
	for(i in 1:length(ids))
		Xr[, i] <- rowSums(X <= X[, ids[i]])
	Xr
}

