setup.mat.samp<-function(cl,type.mt,B=1000,mat.samp=NULL,B.more=0.1,B.max=50000,rand=NA){
	require(multtest)
	if(!is.na(rand))
		set.seed(rand)
	n.cl<-length(cl)
	if(!is.null(mat.samp)){
		if(type.mt=="pairt")
			mat.samp<-pairt.samp.transform(mat.samp)
		if(ncol(mat.samp)!=n.cl)
			stop("The number of columns of mat.samp must be equal to the length of cl.")
		if(type.mt=="f")
			mat.samp<-mat.samp-1		
		a.out<-apply(mat.samp,1,function(a,y=cl) all(sort(a)==sort(y)))
		if(sum(a.out)!=nrow(mat.samp))
			stop("There is something wrong with mat.samp.")
		return(mat.samp)
	}
	if(B<0 || round(B)!=B)
		stop("B must be a positive integer.")
	if(B.more<0 || B.more>1)
		stop("B.more must be between 0 and 1.")
	if(B.max<0 || round(B.max)!=B.max)
		stop("B.max must be a positive integer.")
	B.full<-round(ifelse(type.mt=="pairt",2^(n.cl/2),exp(lgamma(n.cl+1)-sum(lgamma(table(cl)+1)))))
	if(B==0)
		B<-B.full
	B.large<-ceiling(B*(1+B.more))
	if(B.large>=B.full){
		mat.samp<-mt.sample.label(cl,test=type.mt,B=B.large+1)
		cat("\n")
		return(mat.samp)
	}
	if(B.max>=B.full){
		if(type.mt=="pairt"){
			tmp<-pairt.samp.transform(pairt.samp(n.cl/2)) 
	 		cat("We're doing",2^(n.cl/2),"complete permutations \n")
		}	
		else 
			tmp<-mt.sample.label(cl,test=type.mt,B=B.full+1)
		mat.samp<-tmp[sample(B.full,B),]
		cat("and randomly select",B,"of them.\n\n")
		return(mat.samp)
	}
	#if(balanced){
	#	if(!type.mt%in%c("t","t.equalvar"))
	#		stop("The balanced option is only available for the two class unpaired cases.")
	#	n.x<-sum(cl==0)
	#	if(n.x!=n.cl/2 || n.x/2!=round(n.x/2))
	#		stop("Balanced permutations are only possible if n1=n2 are even integer.")
	#	}
	if(type.mt=="pairt"){
		tmp<-matrix(sample(c(-1,1),B*n.cl/2,TRUE),B,n.cl/2)
		mat.samp<-pairt.samp.transform(tmp)
	}
	else{
		mat.samp<-matrix(0,B,n.cl)
		for(i in 1:B) mat.samp[i,]<-sample(cl)
	}
	return(mat.samp)
}
	 
 

