# Copyright (C) 2002 Holger Schwender, University of Dortmund, Germany

# This program calculates the r- and s-values in the computation of the observed test statistics d
# and of the expected d, respectively. This statistics are used in the empBayes-Analysis following
# Efron et al. and in the SAM-Analysis.

# This function could handle missing values and variance zero genes. But we recommend to use some previous
# analysis for a maybe better handling of NAs and variance zero genes. NAs are replaced by the mean of the
# gene. The d-value of genes with variance zero is set to NA.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# data: the used data set 
# x: the columns of the data set which belong to the cases (in the unpaired case) or the "after treatment"- 
#    values (in the paired case) 
# y: the columns of the data set which belong to the control group (in the unpaired case) or to the "before
#    treatment"-measurements
# paired: paired or unpaired data
# mat.samp: The permutation matrix. If specified and correct, this matrix will be used, even if rand and B are
#           specified.
# B: number of permutations used in the calculation of the Null. Will not be used, if mat.samp is specified.
# bal: if TRUE, balanced permutations will be used.
# na.rm: if na.rm=FALSE, the d-values of genes with one or more missing values will be set to NA. If na.rm=TRUE, the
#        missing values will be removed during the computation of the d-values.
# rand: specifies the set.seed for the calculation of the permutation matrix. Will not be used, if mat.samp
#       is specified.


rs.cal<-function(data,x,y=NULL,paired=FALSE,mat.samp=NULL,B=100,bal=FALSE,na.rm=FALSE,rand=NA){
    X<-as.matrix(data[,c(x,y)])  # some adjustments for easier calculations
    mode(X)<-"numeric"
    paired<-ifelse(is.null(y),TRUE,paired)   # new
    NA.genes<-NULL
    if(any(is.na(X))){   # checks if there are NAs
        NA.genes<-unique(ceiling(which(is.na(t(X)))/ncol(X)))   
        cat("Warning: There are",length(NA.genes),"genes with at least one missing value.")
        if(na.rm)
            X[NA.genes,]<-na.replace(X[NA.genes,]) # replace missing values with the gene mean
        if(!na.rm)
            cat(" The d-value of these genes is set to NA.")
        cat("\n","\n")
    }
    if(!paired){      # unpaired case
        n.x<-length(x)    
        n.y<-length(y)
        n<-n.x+n.y
        x.mat<-rep(1,n.x)
        y.mat<-rep(1,n.y)
        mean.x<-as.vector(1/n.x*X[,1:n.x]%*%x.mat)        # computation of the r-values (nominator of d)
        mean.y<-as.vector(1/n.y*X[,(n.x+1):n]%*%y.mat)
        r<-mean.x-mean.y  # r=mean(x)-mean(y)
        center.x<-(X[,1:n.x]-mean.x)^2          # calculation of the s-values (part of the denumerator of d)
        center.y<-(X[,(n.x+1):n]-mean.y)^2
        s<-as.vector(sqrt(((1/n.x+1/n.y)/(n-2))*(center.x%*%x.mat+center.y%*%y.mat)))  
        if(is.null(mat.samp))  # checks if mat.samp is specified
            mat.samp<-sam.sampler(n.x,n.y,B,paired=FALSE,rand=rand,balanced=bal,file.out=NA)   # get the permutation matrix
        if(ncol(mat.samp)!=n)  # checks if mat.samp has the correct number of columns. If not, the function is stopped
            stop("mat.samp has not the correct number of columns.")
        if(!all(mat.samp==1 | mat.samp==0)) # checks if the values of mat.samp are correct. If not, function is stopped
            stop("The values of mat.samp must be 0 or 1.")
        B<-nrow(mat.samp)
        r.perm<-matrix(0,length(r),B)
        s.perm<-matrix(0,length(r),B)
        for(i in 1:B){
            perm<-which(mat.samp[i,]==1)       
            n.perm<-length(perm)      # necessary if balanced permutations are used
            mean.perm.x<-as.vector(1/n.perm*X[,perm]%*%x.mat)      # computation of the r-values for the i-th permutation
            mean.perm.y<-as.vector(1/(n-n.perm)*X[,-perm]%*%y.mat)   # to get the Null    
            r.perm[,i]<-mean.perm.x-mean.perm.y
            center.perm.x<-(X[,perm]-mean.perm.x)^2    # calculation of the s-values for the i-th permutation
            center.perm.y<-(X[,-perm]-mean.perm.y)^2
            s.perm[,i]<-as.vector(sqrt(((1/n.perm+1/(n-n.perm))/(n-2))*(center.perm.x%*%rep(1,n.perm)+center.perm.y%*%rep(1,n-n.perm))))
        Z<-NULL
	}}
    if(paired){  #paired case
        if(!is.null(y) & length(x)!=length(y))  # new # x[i] and y[i] are paired observations. So x and y must have the same length.
            stop("x and y must have the same length.")
        n<-length(x)
        x.mat<-rep(1,n)
        Z<-if(!is.null(y)) X[,1:n]-X[,(n+1):(2*n)]  else X # new # calculation of the r-values of the observed d  (r=mean(x-y))
        r<-as.vector(1/n*Z%*%x.mat)
        center<-(Z-r)^2
        s<-as.vector(sqrt(1/(n*(n-1))*center%*%x.mat))     # calculation of the s-values of the observed d
        if(is.null(mat.samp))       # the same checkings as in the unpaired case
            mat.samp<-sam.sampler(n,n,B,paired=TRUE,rand=rand,balanced=bal,file.out=NA)   #get the permutation matrix
        if(ncol(mat.samp)!=n)
            stop("mat.samp has not the correct number of columns.")
        if(!all(abs(mat.samp)==1))
            stop("The values of mat.samp must be -1 or 1.")
        B<-nrow(mat.samp)
        r.perm<-matrix(0,length(r),B)
        s.perm<-matrix(0,length(r),B)
        for(i in 1:B){
            Z.perm<-t(t(Z)*mat.samp[i,])              # and the same calculations as before for the i-th permutation
            r.perm[,i]<-as.vector(1/n*Z.perm%*%x.mat)
            center.perm<-(Z.perm-r.perm[,i])^2
            s.perm[,i]<-as.vector(sqrt(1/(n*(n-1))*center.perm%*%x.mat))
        }}
    var.0.genes<-NULL
    if(any(s==0,na.rm=TRUE)){    # sets the values of r and s of each gene with variance 0 to NA
        cat("Warning: There are",sum(s==0,na.rm=TRUE),"genes which have variance Zero or no non-missing values.","\n",
            "        The d-value of these genes is set to NA.","\n","\n")
        var.0.genes<-which(s==0)
        r[var.0.genes]<-NA
        s[var.0.genes]<-NA
        r.perm[var.0.genes,]<-NA
        s.perm[var.0.genes,]<-NA
    }

    structure(list(r=r,s=s,r.perm=r.perm,s.perm=s.perm,Z=Z,mat.samp=mat.samp,var.0.genes=var.0.genes,
		NA.genes=NA.genes))
}
