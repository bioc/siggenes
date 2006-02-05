rs.cal<-function(data,x,y=NULL,paired=FALSE,mat.samp=NULL,B=100,bal=FALSE,na.rm=FALSE,rand=NA){
    	X<-as.matrix(data[,c(x,y)]) 
    	mode(X)<-"numeric"
    	paired<-ifelse(is.null(y),TRUE,paired)   
    	NA.genes<-NULL
    	if(any(is.na(X))){   
        	NA.genes<-unique(ceiling(which(is.na(t(X)))/ncol(X)))   
        	cat("Warning: There are",length(NA.genes),"genes with at least one missing value.")
        	if(na.rm)
            		X[NA.genes,]<-na.replace.cont(X[NA.genes,])
        	if(!na.rm)
            		cat(" The d-value of these genes is set to NA.")
        	cat("\n","\n")
    	}
    	if(!paired){
        	n.x<-length(x)    
        	n.y<-length(y)
        	n<-n.x+n.y
        	mean.x<-rowMeans(X[,1:n.x])
        	mean.y<-rowMeans(X[,(n.x+1):n])
        	r<-mean.x-mean.y  
        	center.x<-(X[,1:n.x]-mean.x)^2  
        	center.y<-(X[,(n.x+1):n]-mean.y)^2
        	s<-sqrt(((1/n.x+1/n.y)/(n-2))*(rowSums(center.x)+rowSums(center.y)))  
        	if(is.null(mat.samp))  
            		mat.samp<-sam.sampler(n.x,n.y,B,paired=FALSE,rand=rand,balanced=bal,
				file.out=NA)
        	if(ncol(mat.samp)!=n) 
            		stop("mat.samp has not the correct number of columns.")
        	if(!all(mat.samp==1 | mat.samp==0)) 
            		stop("The values of mat.samp must be 0 or 1.")
        	B<-nrow(mat.samp)
        	r.perm<-matrix(0,length(r),B)
        	s.perm<-matrix(0,length(r),B)
        	for(i in 1:B){
            		perm<-which(mat.samp[i,]==1)       
            		n.perm<-length(perm)
            		mean.perm.x<-rowMeans(X[,perm]) 
            		mean.perm.y<-rowMeans(X[,-perm])    
            		r.perm[,i]<-mean.perm.x-mean.perm.y
            		center.perm.x<-(X[,perm]-mean.perm.x)^2 
            		center.perm.y<-(X[,-perm]-mean.perm.y)^2
            		s.perm[,i]<-sqrt(((1/n.perm+1/(n-n.perm))/(n-2))
				*(rowSums(center.perm.x)+rowSums(center.perm.y)))
        		Z<-NULL
		}
	}	
    	if(paired){ 
        	if(!is.null(y) & length(x)!=length(y)) 
            		stop("x and y must have the same length.")
        	n<-length(x)
        	x.mat<-rep(1,n)
        	Z<-if(!is.null(y)) X[,1:n]-X[,(n+1):(2*n)]  else X 
        	r<-rowMeans(Z)
        	center<-(Z-r)^2
        	s<-sqrt(1/(n*(n-1))*rowSums(center)) 
	        if(is.null(mat.samp))    
		mat.samp<-sam.sampler(n,n,B,paired=TRUE,rand=rand,balanced=bal,file.out=NA)
        	if(ncol(mat.samp)!=n)
            		stop("mat.samp has not the correct number of columns.")
        	if(!all(abs(mat.samp)==1))
            		stop("The values of mat.samp must be -1 or 1.")
        	B<-nrow(mat.samp)
       	 	r.perm<-matrix(0,length(r),B)
        	s.perm<-matrix(0,length(r),B)
        	for(i in 1:B){
           		Z.perm<-t(t(Z)*mat.samp[i,])  
            		r.perm[,i]<-rowMeans(Z.perm)
            		center.perm<-(Z.perm-r.perm[,i])^2
            		s.perm[,i]<-sqrt(1/(n*(n-1))*rowSums(center.perm%*%x.mat))
        	}
	}
    	var.0.genes<-NULL
    	if(any(s==0,na.rm=TRUE)){    
        	cat("Warning: There are",sum(s==0,na.rm=TRUE),
			"genes which have variance Zero or no non-missing values.","\n",
            		"        The d-value of these genes is set to NA.","\n","\n")
        	var.0.genes<-which(s==0)
        	r[var.0.genes]<-NA
        	s[var.0.genes]<-NA
        	r.perm[var.0.genes,]<-NA
        	s.perm[var.0.genes,]<-NA
    	}
	structure(list(r=r,s=s,r.perm=r.perm,s.perm=s.perm,Z=Z,mat.samp=mat.samp,
		var.0.genes=var.0.genes,NA.genes=NA.genes))
}
 
 

