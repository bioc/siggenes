# Copyright (c) 2002 Holger Schwender

# This function generates the so called Roller Coaster Plots, i.e. the plot of delta vs. FDR and the plot
# of delta vs. #significant genes.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# mat.fdr: a matrix with the interesting statistics (i.e. delta, #number of significant genes and FDR)
#          which is, e.g., provided by sam.fdr()
# helplines: if TRUE, helplines are produced in the plots for an easier evaluation



roller.coaster<-function(mat.fdr,helplines=TRUE){
    par(mfrow=c(1,2),lab=c(10,8,7))  # two plots on one graphsheet
    invisible() 
    plot(mat.fdr[,"delta"],100*mat.fdr[,"fdr"],main="Delta vs. FDR",xlab="delta",ylab="FDR (in %)",
        type="b")   # delta vs.FDR
    if(helplines){  # helplines are plotted
        segments(0,100*mat.fdr[,"fdr"],mat.fdr[,"delta"],100*mat.fdr[,"fdr"],lty=2)
        segments(mat.fdr[,"delta"],100*mat.fdr[,"fdr"],mat.fdr[,"delta"],-100*max(mat.fdr[,"fdr"]),lty=2)
    }
    plot(mat.fdr[,"delta"],mat.fdr[,"nsig"],main="Delta vs. Significant Genes",xlab="delta",
        ylab="number of significant genes",type="b")  # delta vs. #significant genes
    if(helplines){ # helplines are plotted
        segments(0,mat.fdr[,"nsig"],mat.fdr[,"delta"],mat.fdr[,"nsig"],lty=2)
        segments(mat.fdr[,"delta"],mat.fdr[,"nsig"],mat.fdr[,"delta"],-max(mat.fdr[,"nsig"]),lty=2)
    }
    par(mfrow=c(1,1))
}
