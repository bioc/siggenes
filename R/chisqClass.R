"chisqClass" <-
function(data,cl,n.cat,check=TRUE,pam=FALSE){
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(missing(n.cat))
		n.cat<-max(data)
	uni.cl<-sort(unique(cl))
	n.lev<-length(uni.cl)
	if(length(n.lev)>10)
		stop("cl contains more than 10 different values.")
	if(any(uni.cl!=1:n.lev))
		stop("The labels of the classes must be 1 to ",n.lev,".")
	n.obs<-ncol(data)
	n.snp<-nrow(data)
	if(length(cl)!=n.obs)
		stop("The length of cl must be equal to the number of observations.")
	CL<-matrix(0,n.obs,n.lev)
	for(i in 1:n.lev)
		CL[cl==i,i]<-1
	vec.ncl<-colSums(CL)
	listCells<-vector("list",n.cat)
	for(i in 1:n.cat)
		listCells[[i]]<-data==i
	if(check){
		for(i in 1:n.cat){
			tmp.rS<-rowSums(listCells[[i]])
			if(any(tmp.rS==0))
				stop("All variables must have the same number of levels.")
		}
	}
	if(pam){
		mat.stats<-matrix(0,n.snp,n.lev)
		list.Nobs<-vector("list",n.cat)
		list.Nexp<-vector("list",n.cat)
		for(i in 1:n.cat){
			tmp<-computeContCols(listCells[[i]],CL,vec.ncl,n.obs,pam=TRUE)
			mat.stats<-mat.stats+tmp$mat.stat
			list.Nobs[[i]]<-tmp$N.obs
			list.Nexp[[i]]<-tmp$N.exp
		}
		tmp<-t(vec.ncl)%x%rep(1,n.snp)
		out<-mat.stats-tmp
		if(!is.null(rownames(data)))
			rownames(out)<-rownames(data)
		return(list(mat.stat=out,list.Nobs=list.Nobs,list.Nexp=list.Nexp,n=n.obs))
	}
	else{
		listNexp<-lapply(listCells,computeContCols,CL=CL,vec.ncl=vec.ncl,n.obs=n.obs)
		mat.stats<-matrix(unlist(listNexp),ncol=n.cat)
		out<-rowSums(mat.stats)-n.obs
	}
	if(!is.null(rownames(data)))
		names(out)<-rownames(data)
	out
}

