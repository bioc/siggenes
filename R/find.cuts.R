# Copyright (c) 2002 Holger Schwender

# This program finds the cutup and the cutlow in the SAM Analysis. It also calculates the FDR.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# delta: for this value of delta the cutup and cutlow will be computed
# d.sort: the vector of the sorted observed d-values
# d.diff: the vector of the difference of the d-values and the d.bar-values where the d.bar-vector contains
#         the expected d-values and the d-vector consists of the observed d-values
# d.perm: matrix of permuted d-values. The rowwise average of d.perm is d.bar. In the sam.wilc() case d.perm
#         is a matrix with 2 columns with the possible values for W in the first and the corresponing 
#         expected number of W-values in the second column
# p0: probability that a gene is unaffected
# j0: the index of the d.bar value which is closest to 0
# med: if med=TRUE, the median number of falsely called genes will be used in the calculation of the FDR. 
#      Otherwise the expected number will be used
# wilc: if TRUE, the computation will be done for sam.wilc(). If FALSE, it will be done for sam().

find.cuts<-function(delta,d.sort,d.diff,d.perm,p0,j0,med,wilc=FALSE){
    m<-length(na.exclude(d.sort))  # number of genes
    # the index j1 > j0 of the gene is computed which is the first gene starting from the origin and going
    # to the right for which d.diff > delta  
    j1<-ifelse(any(d.diff[j0:m]>=delta),j0-1+min(which(d.diff[j0:m]>=delta)),m+1)
    # cutup is calculated; if there is no j1 so that d.diff[j1] > delta, cutup will be set to Inf
    cutup<-ifelse(j1!=m+1,d.sort[j1],Inf) 
    # the index j2 < j0 of the gene is computed which is the first gene starting from the origin and going
    # to the left for which d.diff < -delta 
    j2<-ifelse(any(d.diff[1:j0]<= -delta),max(which(d.diff[1:j0]<= -delta)),0)
    # cutlow is calculated; if there is no j2 so that d.diff[j2] < -delta, cutlow will be set to -Inf
    cutlow<-ifelse(j2!=0,d.sort[j2],-Inf)
    nsig<-m-j1+1+j2 # number of significant genes
    if(!wilc){
        if(!med)
            false<-sum(d.perm>=cutup | d.perm<=cutlow,na.rm=TRUE)/ncol(d.perm)  # expected number of falsely called genes
        if(med){
            vec.false<-NULL
            for(i in 1:ncol(d.perm))
                vec.false[i]<-sum(d.perm[,i]>=cutup | d.perm[,i]<=cutlow,na.rm=TRUE)
            false<-median(vec.false)   # median number of falsely called genes
    }}
    if(wilc)
        false<-sum(d.perm[which(d.perm[,1]>=cutup | d.perm[,1]<=cutlow),2])
    fdr<-p0*false/max(nsig,1)    # FDR
    cut.out<-c(delta=delta,p0=p0,false=false,nsig=nsig,fdr=fdr,cutlow=cutlow,cutup=cutup,j1=j1,j2=j2)
    invisible(cut.out)
}
