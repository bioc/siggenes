# Copyright (c) 2002 Holger Schwender

# This function computes for a given number or proportion of genes a delta so that about that number or
# proportion of genes fall outside this thresholds cutlow and cutup, i.e. are called significant. 
# Sometimes it is not possible that exactly this number or proportion of genes fall outside delta. If such
# a situation occurs, a lower and upper bound will be given.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# d.sort: vector of sorted observed d-values
# d.bar: vector of sorted expected d-values
# d.perm: matrix of sorted permuted d-values. The rowwise mean of d.perm is d.bar
# ngenes: number or percentage of genes which should fall outside delta, i.e. which should be called significant
# iteration: number of iterations which will be used in the search for delta. Default is 3. This should usually
#            work.
# initial.delta: a vector of initial guesses for delta. Default is (.2, .4, .6, ..., 1.8, 2)
# med: if TRUE, the median number of falsely called genes will be computed. If FALSE, the expected number will
#      be calculated
# p0: the (prior) probability that a gene is unaffected


sam.ngenes<-function(d.sort,d.bar,d.perm,ngenes=0.05,iteration=3,initial.delta=seq(.2,2,.2),med=TRUE,p0=NA){
    if(ngenes!=round(ngenes)){
        if(ngenes<=0 || ngenes>=1)
            stop("ngenes must be an positive integer or between 0 and 1")
        ngenes<-floor(ngenes*length(d.sort))  # the proportion of genes is transformed to the number of genes
    }
    if(ngenes<=0)
        stop("ngenes must be an positive integer")
    cat("SAM Analysis for",ngenes,"genes:","\n")
    delta.ngenes<-initial.delta
    for(i in 1:iteration){     # for a vector of delta some statistics (e.g. #significant genes) is computed
        fdr.ngenes<-sam.fdr(d.sort,d.bar,d.perm,p0,delta=delta.ngenes,med=med,graphic=FALSE)$mat.fdr
        if(any(fdr.ngenes[,4]==ngenes)){    # if a delta is found which leads to ngenes significant genes
            fdr.ngenes<-as.vector(fdr.ngenes[min(which(fdr.ngenes[,4]==ngenes)),1:5])   # the search is done
            delta.ngenes<-fdr.ngenes[1]                                      # and some output is made
            names(fdr.ngenes)<-c("delta","p0","false","called","FDR")
            cat("Set delta =",round(delta.ngenes,4),"to get",ngenes,"significant genes.","\n","\n")
            print(fdr.ngenes)
            invisible(return(list(fdr.ngenes=fdr.ngenes,delta.ngenes=delta.ngenes,p0=p0)))
        }
        # if no such delta was found, a new vector of deltas is been computed. The minimum delta is that
        # delta which leads to the nearest smaller number of significant genes of the previous vector of delta.
        # The max delta is that delta that leads to the nearest higher number of significant genes. 
        delta.ngenes<-seq(fdr.ngenes[max(which(fdr.ngenes[,4]>ngenes)),1],fdr.ngenes[min(which(fdr.ngenes[,4]<ngenes)),1],
            length=20)
    }
    # Sometimes it is not possible to find such a delta. Then a lower and upper bound are given.
    cat("It is not possible to determine a delta for exactly",ngenes,"genes.","\n","\n",
        "Lower and upper bound:","\n")
    fdr.ngenes<-as.data.frame(rbind(lower=fdr.ngenes[max(which(fdr.ngenes[,4]>ngenes)),],
        upper=fdr.ngenes[min(which(fdr.ngenes[,4]<ngenes)),]))
    print(fdr.ngenes)
    delta.ngenes<-NULL
    invisible(structure(list(fdr.ngenes=fdr.ngenes,delta.ngenes=delta.ngenes,p0=p0)))
}
    
    
    
    
