# Copyright (c) 2002 Holger Schwender

# This function computes the observed and the expected Wilcoxon Rank statistics and the distribution of
# the W-values under the Null.

# This function is a helpfunction for sam.wilc() and ebam.wilc().

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# data: the used data set. Every row of this data set must correspond to a gene.
# x: the columns of data, which belong to the cases (unpaired) or to the "after treatment"-measurements (paired)
# y: the columns of data, which belong to the control group (unpaired) or to the "before treatment"-measurements 
#    (paired)
# paired: paired or unpaired data. If paired=TRUE, i.e. there are paired observations, x and y must have the same
#         length and (x[i],y[i]) must be an observation pair
# zero.rand: if there are Zeros in the paired case, what should be done with these Zeros? If zero.rand=TRUE, a sign
#            is randomly assigned to each observation pair (x[i],y[i]) with x[i]-y[i]=0. If FALSE, the sign of such
#            observation pairs is set to '-' (suggested by Lehmann (1975))
# rand: is only necessary in the paired case for the random choice of the signs of Zeros (see zero.rand)
# na.rm: if na.rm=TRUE, each missing value in the data will be replaced by the genewise mean. If FALSE, the W-value
#        of every gene with one or more missing values is set to NA


wilc.cal<-function(data,x,y,paired=FALSE,zero.rand=TRUE,rand=NA,na.rm=FALSE){
    if(!is.na(rand))
        set.seed(rand)
    X<-as.matrix(data[,c(x,y)])
    mode(X)<-"numeric"
    n.genes<-nrow(X)  # number of genes
    NA.genes<-NULL       
    var.0.genes<-NULL
    if(any(is.na(X))){  # checks if there are NAs
        NA.genes<-unique(ceiling(which(is.na(t(X)))/ncol(X)))  # which gene has NAs?
        cat("Warning: There are",length(NA.genes),"genes with at least one missing value.")
        if(na.rm){     
            X[NA.genes,]<-na.replace(X[NA.genes,])   # replace missing values with the gene mean
            X[which(is.na(X))]<-0  # if there are still NAs, i.e. if there are genes with no non-missing
        }                          # value, these NAs will be set to 0 and treated later.
        if(!na.rm)
            cat(" The W-value of these genes is set to NA.")
        cat("\n","\n")
    }
    
    
    if(!paired){  # paired case
        n.x<-length(x)   
        n.y<-length(y)
        n<-n.x+n.y
        W.mean<-n.x*(n+1)/2   # mean of W-values under the null
        W.min<-n.x*(n.x+1)/2  # minimum of W-values under the null
        W.max<-n.x*(2*n.y+n.x+1)/2  # max of W-values under the null
        X.rank<-t(apply(X,1,rank))  # compute the rowwise ranks
        X.sum<-apply(X.rank,1,var)  # check if there are some genes with variance Zero
        if(any(X.sum==0)){     # which genes have variance zero?
            cat("Warning: There are",sum(X.sum==0),"genes with variance Zero or no non-missing value. The W-value of these genes is set to NA.",
                "\n","\n")
            var.0.genes<-which(X.sum==0)        
        }
        W<-rowSums(X.rank[,1:n.x])  # compute the observed W-value of each gene
        if(!is.null(var.0.genes)){
            W[var.0.genes]<-NA            # set the W-value of the genes with variance zero to NA
            n.genes<-n.genes-length(var.0.genes)
        }
        if(!na.rm && !is.null(NA.genes)){
            W[NA.genes]<-NA             # set the W-value of the genes which have variance zero to NA
            n.genes<-n.genes-length(NA.genes)   # number of genes with non-missing W-value
        }
        W.exp<-W.min+qwilcox(((1:n.genes)-.5)/n.genes,n.x,n.y)  # compute the expected W-values
        f.null<-dwilcox(0:(n.x*n.y),n.x,n.y)   # compute the distribution of the W-values under the Null
    }
    if(paired){  # paired case
        if(length(x)!=length(y))  # check if x and y have the same length
            stop("x any y must have the same length.")
        n<-length(x)
        X<-X[,1:n]-X[,(n+1):(2*n)]   # subtract the y-columns from the corresponding x-columns of the data set 
        X.sum<-apply(X,1,var)     # check if there are genes with variance zero
        if(any(X.sum==0)){   
            cat("There are",sum(X.sum==0,na.rm=TRUE),"genes with variance Zero or no non-missing value. The W-value of these genes is set to NA.",
                "\n","\n")
            var.0.genes<-which(X.sum==0)  # which genes have variance zero?
        }
        W.max<-n*(n+1)/2   # max of W-values under the Null, min is always 0
        W.mean<-n*(n+1)/4  # mean of W-values under the Null
        if(sum(X==0)>0){   # check if there are some Zeros
            cat("There are",sum(X==0),"Zeros.","\n","\n")
            if(zero.rand)   # if zero.rand=TRUE, a sign is randomly assigned to observation pairs with Zero difference 
                X[which(X==0)]<-sample(c(.00001,-.00001),sum(X==0),rep=TRUE)
        }
        W<-NULL
        X.rank<-NULL
        for(i in 1:n.genes)   # compute the observed W-values
            W[i]<-sum(rank(abs(X[i,]))[X[i,]>0])
        if(!is.null(var.0.genes)){  
            W[var.0.genes]<-NA    # set the W-values of genes with variance Zero to NA       
            n.genes<-n.genes-length(var.0.genes)
        }
        if(!na.rm && !is.null(NA.genes)){
            W[NA.genes]<-NA      # set the W-values of genes with missing values to NA
            n.genes<-n.genes-length(NA.genes)   # number of genes with non-missing W-values
        }
        W.exp<-qsignrank(((1:n.genes)-.5)/n.genes,n) # expected W-values
        f.null<-dsignrank(0:W.max,n)    # compute the distribtuion of the W-values under the null   
    }
    if(sum(W!=round(W),na.rm=TRUE)>0){  
        cat("tied Wilcoxon scores:", sum(W!=round(W),na.rm=TRUE),"\n","\n")
        y.rand<-sample(c(-0.5,0.5),length(which(W!=round(W))),replace=TRUE)            # integer
        W[which(W!=round(W))]<-W[which(W!=round(W))]+y.rand
    }
    names(f.null)<-as.character(ifelse(paired,0,W.min):W.max) 
    
    structure(list(W=W,W.exp=W.exp,f.null=f.null,X.rank=X.rank,var.0.genes=var.0.genes,
		NA.genes=NA.genes,n.genes=n.genes))
}
    
