# Copyright (c) 2002 Holger Schwender

# This function computes the p-values for the Wilcoxon (Sign) Rank Statistics which were calculated in 
# sam.wilc().

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# W: the observed unsorted W-values
# p0: the probability that a gene is unaffected
# n.x: length of x (see sam.wilc())
# n.y: length of y (see sam.wilc())
# paired: paired or unpaired data

q.value.wilc<-function(W,p0,n.x,n.y,paired=FALSE){
    m<-length(na.exclude(W))   # number of genes with non missing W-value  
    n<-ifelse(paired,n.x,n.x+n.y)  # number of independent observations
    W.mean<-ifelse(paired,n*(n+1)/4,n.x*(n+1)/2)   # mean of W-values under the null
    W.min<-ifelse(paired,0,n.x*(n.x+1)/2)  # minimum of W-values under the null
    W.max<-ifelse(paired,n*(n+1)/2,n.x*(2*n.y+n.x+1)/2)  # max of W-values under the null
    p.value<-numeric(W.max-W.min+1)  # vector of p-values
    nsig<-numeric(W.max-W.min+1)     # vector of significant genes
    for(i in W.min:floor(W.mean)){  # calculation of the p-values and the corresponding number of significant genes
        nsig[i-W.min+1]<-sum(W<=i | W>=W.max+W.min-i,na.rm=TRUE)
        nsig[length(nsig)+W.min-i]<-sum(W<=i | W>=W.max+W.min-i,na.rm=TRUE)
        p.value[i-W.min+1]<-min(1,2*ifelse(paired,psignrank(i,n),pwilcox(i-W.min,n.x,n.y)))
        p.value[length(p.value)+W.min-i]<-min(1,2*ifelse(paired,psignrank(i,n),pwilcox(i-W.min,n.x,n.y)))
    }
        
    numeric(W.max-W.min+1)->vec.q.value      # vector of possible q-values
    vec.q.value[floor(W.mean-W.min+1)]<-p0*p.value[floor(W.mean-W.min+1)]  # initial setting for the q-values
    for(i in (floor(W.mean)-1):W.min)   # computation of the q-values
        vec.q.value[i-W.min+1]<-min(p0*p.value[i-W.min+1]/(nsig[i-W.min+1]/m),vec.q.value[i-W.min+2])
    vec.q.value[(ceiling(W.mean):W.max)-W.min+1]<-rev(vec.q.value[(W.min:floor(W.mean))-W.min+1])
    q.value<-numeric(length(W))  # vector of gene-specific q-values 
    for(i in W.min:W.max)
        q.value[which(W==i)]<-vec.q.value[i-W.min+1]
    mat.qvalue<-cbind(W,q.value)[order(abs(W-W.mean),na.last=TRUE),]   # matrix of ordered W- and corresponding q-values
    structure(list(mat.qvalue=mat.qvalue,p.value=p.value,nsig=nsig,q.value=q.value,
		vec.q.value=vec.q.value,W.min=W.min,W.max=W.max,W.mean=W.mean))
}
