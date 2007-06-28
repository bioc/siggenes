recodeVal <- function(mat, lev){
	if(!is.list(lev))
		stop("lev must be a list.")
	n.lev <- unlist(lapply(lev, length))
	n.lev <- unique(n.lev)
	if(length(n.lev) != 1)
		stop("Each vector in lev must have the same length.")
	for(i in 1:length(lev)){
		tmp <- lev[[i]]
		mat.lev <- matrix(FALSE, nrow(mat), n.lev)
		for(j in 1:n.lev)
			mat.lev[,j] <- rowSums(mat == tmp[j])
		ids <- rowSums(mat.lev) == ncol(mat)
		tmp.mat <- mat[ids, ]
		for(j in 1:n.lev)
			tmp.mat[tmp.mat == tmp[j]] <- j
		mat[ids, ] <- tmp.mat
	}
	mat
}


recodeLevel<-function(mat,lev){
	if(is.list(lev))
		return(recodeVal(data,lev))
	for(i in 1:length(lev))
		data[data==lev[i]]<-i
	data
}


		 