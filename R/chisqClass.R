"chisqClass" <-
function(data,cl,n.cat,check=TRUE, withNA=FALSE){
	#if(any(is.na(data)))
	#	stop("No NAs allowed.")
	if(missing(n.cat))
		n.cat<-max(data, na.rm=TRUE)
	uni.cl<-sort(unique(cl))
	n.lev<-length(uni.cl)
	if(length(n.lev)>10)
		stop("cl contains more than 10 different values.")
	if(any(uni.cl!=1:n.lev))
		stop("The labels of the classes must be 1 to ",n.lev,".")
	#n.obs <- if(!withNA) ncol(data) else rowSums(!is.na(data))
	#n.snp<-nrow(data)
	#if(length(cl)!=n.obs)
	#	stop("The length of cl must be equal to the number of observations.")
	CL <- matrix(0, ncol(data) , n.lev)
	for(i in 1:n.lev)
		CL[cl==i,i]<-1
	listCells<-vector("list",n.cat)
	if(!withNA){
		vec.ncl<-colSums(CL)
		n.obs <- ncol(data)
		for(i in 1:n.cat)
			listCells[[i]] <- data==i
	}
	else{
		mat.notna <- !is.na(data)
		vec.ncl <- mat.notna%*%CL
		n.obs <- rowSums(mat.notna)
		for(i in 1:n.cat)
			listCells[[i]] <- mat.notna & data==i
	}
	if(check){
		for(i in 1:n.cat){
			tmp.rS<-rowSums(listCells[[i]])
			if(any(tmp.rS==0))
				stop("All variables must show the same number of levels.")
		}
	}
	listNexp<-lapply(listCells,computeContCols,CL=CL,vec.ncl=vec.ncl,n.obs=n.obs)
	mat.stats<-matrix(unlist(listNexp),ncol=n.cat)
	out<-rowSums(mat.stats)-n.obs
	if(!is.null(rownames(data)))
		names(out)<-rownames(data)
	out
}

