na.handling<-function(X,na.replace=TRUE,na.method="mean"){
	if(!is.matrix(X))
		stop("X must be a matrix.")
	vec.na<-rowSums(is.na(X))
	NA.genes<-which(vec.na>0)
	if(length(NA.genes)==0)
		return(list(X=X,NA.genes=NULL))
	n.col<-ncol(X)
	if(!na.replace){
		X<-X[-NA.genes,]
		warning("There are ",length(NA.genes)," genes with at least one missing expression value.",
			"\n","All these genes are removed, and their d-values are set to NA.",call.=FALSE)
	}
	else{
		X[NA.genes,]<-na.replace.cont(X[NA.genes,],na.method=na.method)
		warning("There are ",length(NA.genes)," genes with at least one missing expression value.",
			"\n","The NAs are replaced by the gene-wise ",na.method,".",call.=FALSE)
		if(any(vec.na>=n.col-1)){
			NA.genes0<-which(vec.na==n.col)
			NA.genes1<-which(vec.na==n.col-1)
			warning(length(NA.genes0)," of the ",length(NA.genes)," genes with at least one NA have no and ",
				length(NA.genes1)," have one non-missing expression value.","\n","All these ",
				length(c(NA.genes0,NA.genes1))," genes are removed, and their d-values are set to NA.",
				call.=FALSE)
			NA.genes<-c(NA.genes0,NA.genes1)
			X<-X[-c(NA.genes0,NA.genes1),]
		}
		else
			NA.genes<-NULL
	}
	list(X=X,NA.genes=NA.genes)
}
 
 

