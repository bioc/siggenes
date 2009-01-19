rowRanksWilc <- function(X, y){
	if(length(y)!=ncol(X))
		stop("The length of y must be equal to the number of columns of X.")
	ids <- which(y == 1)
	Xr <- matrix(0, nrow(X), length(ids))
	for(i in 1:length(ids))
		Xr[, i] <- rowSums(X <= X[, ids[i]])
	Xr
}

