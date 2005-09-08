cat.stat<-function(data,cl,B=1000,med=FALSE,gene.names=NULL,na.replace=TRUE,check.levels=TRUE,
		rand=NA){
	if(!is.na(rand))
		set.seed(rand)
	n.cl<-length(cl)
	n.snps<-nrow(data)
	uni.cl<-unique(cl)
	msg<-paste("SAM Analysis for Categorical Data with",length(uni.cl),"Classes \n\n")
	if(!is.null(gene.names))
		gene.names<-gene.names
	if(ncol(data)!=n.cl)
		stop("The length of cl must be equal to the number of columns of data.")
	if(length(uni.cl)==1)
		stop("There must be at least two groups.")
	if(any(table(cl)<10))
		stop("There must be at least 10 observations in each group.")
	if(n.cl<25)
		stop("There should be at least 25 samples.")
	if(check.levels){
		n.lev<-numeric(n.snps)
		for(i in 1:n.snps)
			n.lev[i]<-length(table(data[i,]))
		if(any(n.lev>=10))
			stop("Some SNPs have 10 or more levels.")
		if(length(unique(n.lev))!=1)
			stop("All SNPs must have the same number of levels.")
	}
	if(any(rowSums(is.na(data))>0)){
		n.na<-rowSums(is.na(data))
		NA.genes<-which(n.na>0)
		warning("There are ",length(NA.genes)," SNPs with at least one missing value.",
			if(na.replace) " These NAs will be replaced \n",
			" by random draws from the distribution of the corresponding SNP.",call.=FALSE)
		if(any(n.na>n.cl-25)){
			many.na<-which(n.na>n.cl-25)
			warning(length(many.na)," of the ",length(NA.genes)," SNPs with at least one ",
				"NA have less than 25 non-missing values.",
				"\n","All these SNPs are removed.",call.=FALSE)
			data<-data[-many.na,]
			n.na<-n.na[-many.na]
			NA.genes<-NA.genes[!NA.genes%in%many.na]
		}
		if(na.replace & length(NA.genes)!=0)
			data[NA.genes,]<-if(length(NA.genes)==1) na.replace.dist(data[NA.genes,])
				else t(apply(data[NA.genes,],1,na.replace.dist))		
	}
	nr<-length(uni.cl)
	n.row<-nrow(data)
	d<-numeric(n.row)
	d.perm<-matrix(0,n.row,B)
	for(i in 1:n.row){
		if(!is.na(rand))
			set.seed(rand)
		tab.chi<-table(cl,data[i,])
		n.chi<-sum(tab.chi)
		sr<-rowSums(tab.chi)
		sc<-colSums(tab.chi)
		e.chi<-outer(sr,sc)/n.chi
		d[i]<-sum((tab.chi-e.chi)^2/e.chi)
		nc<-ncol(tab.chi)
		d.perm[i,]<-.C("chisqsim",as.integer(nr),as.integer(nc),as.integer(sr),
			as.integer(sc),as.integer(n.chi),as.integer(B),as.double(e.chi),
			integer(nr*nc),double(n.chi+1),integer(nc),results=double(B),
			PACKAGE="stats")$results
	}
	d.perm<-apply(d.perm,2,sort)
	d.bar<-rowMeans(d.perm)
	d.rank<-rank(-d,ties="first")
	mat.rank<-matrix(0,n.row,B)
	for(i in 1:B)
		mat.rank[,i]<-rank(-c(d.perm[,i],d),ties="first")[n.row+(1:n.row)]-d.rank
	vec.false<-if(med) apply(mat.rank,1,median) else rowMeans(mat.rank)
	p.value<-(1+rowSums(mat.rank))/(n.row*B+1)
	if(n.snps!=n.row){
		d.new<-vec.new<-p.new<-rep(NA,n.snps)
		d.new[-many.na]<-d
		vec.new[-many.na]<-vec.false
		p.new[-many.na]<-p.value
		d<-d.new
		vec.false<-vec.new
		p.value<-p.new
	}
	rm(d.perm)
	if(!is.null(gene.names))
		names(d)<-names(p.value)<-names(vec.false)<-substring(gene.names,1,50)
	structure(list(d=d,d.bar=d.bar,p.value=p.value,vec.false=vec.false,discrete=FALSE,
		s=numeric(0),s0=numeric(0),mat.samp=matrix(numeric(0)),msg=msg,fold=numeric(0)))		
	
} 
 

