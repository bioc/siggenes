# Copyright (c) 2002 Holger Schwender

# This function computes the FDR for several delta. The output of this program is a table with some
# statistics (p0, number of significant genes, number of falsely called genes, FDR) for these delta
# and a SAM plot of the expected vs. the observed d-values for some delta (not necessary the same delta
# as in the table). There will also be a plot of delta vs. FDR and delta vs. #significant genes in the
# output. This table and these plots should help to choose delta.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# d.sort: vector of sorted observed d-values
# d.bar: vector of sorted expected d-values
# d.perm: matrix of sorted permuted d-values. The rowwise mean of this matrix is d.bar
# p0: the (prior) probability that a gene is unaffected
# delta: vector of values for delta for which the statistics (see above) should be computed
# med: if TRUE, the median number of falsely called genes will be computed. If FALSE, the expected number will
#      be calculated
# graphic: if TRUE, a SAM plot of the expected vs. the observed d-values for some delta (see thres) will be made
#          as well as a plot of delta vs. FDR and delta vs. #significant genes 
# thres: vector of the delta for which the SAM Plot is made
# pty.square: if TRUE, a square SAM Plot is generated with x- and y-axes having the same range
# helplines: if TRUE, helplines are plotted in both the plot of delta vs. FDR and the plot of
#            delta vs. #significant genes
# wilc: if TRUE, an analysis is done for sam.wilc(). Otherwise (default) for sam().


sam.fdr<-function(d.sort,d.bar,d.perm,p0,delta=(1:10)/5,med=TRUE,graphic=TRUE,thres=seq(.5,2,.5),
        pty.square=TRUE,helplines=TRUE,wilc=FALSE){
    d.diff<-d.sort-d.bar  # calculation of the difference of the observed and the corresponding expected
                               # d-values
    j0<-ifelse(!wilc,which(abs(d.bar)==min(abs(d.bar),na.rm=TRUE)),floor(length(na.exclude(d.bar))/2)) 
             # interpretation of Tushers "start at the origin", the index of
    mat.fdr<-NULL                           # the expected d-value which is closest to 0 is computed
    for(i in 1:length(delta)){
        cuts.out<-find.cuts(delta[i],d.sort,d.diff,d.perm,p0,j0,med,wilc=wilc)    # calculation of cutlow, #significant genes,
        mat.fdr<-rbind(mat.fdr,cuts.out[-1])           # #falsely called genes, FDR
    }
    if(length(delta)==1){  # distinction is necessary for further analysis
        mat.fdr<-matrix(c(delta,mat.fdr),1)
        tab.fdr<-as.vector(round(mat.fdr[,1:5],3))
    }
    if(length(delta)>1){
        mat.fdr<-cbind(delta,mat.fdr)
        tab.fdr<-as.data.frame(round(mat.fdr[,1:5],3))
    }
    names(tab.fdr)<-c("delta","p0","false","called","FDR")
    if(graphic){      # SAM Plot
        sam.plotter(d.sort,d.bar,thres,pty.square=pty.square, main="SAM Plot for some delta",
            color=2:(length(thres)+1),make.legend=TRUE)  # SAM Plot
        X11()   
        roller.coaster(mat.fdr,helplines=helplines) # Delta vs. FDR and vs. #significant
                                # genes, respectively
    }
    structure(list(tab.fdr=tab.fdr,mat.fdr=mat.fdr,p0=p0))
}
