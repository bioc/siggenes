# Copyright (c) 2002 Holger Schwender

# This function makes a SAM analysis using Wilcoxon Rank Statistics. The output of this function is
# a table of statistics (p0, #number of significant genes, #number of falsely called genes, #FDR) for some
# delta, a SAM Plot for some (not necessarily the same) delta and the plots of delta vs. FDR and 
# delta vs. #significant genes.

# The output of this function can be used to select the Delta in the SAM Analysis. For further analysis
# with sam.plot() the output of sam.wilc() must be assigned to an object. 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# data: the data set. Every row of this data set must represent a gene.
# x: the columns of the data set which belong to the cases (in the unpaired case) or the "after treatment"- 
#    values (in the paired case) 
# y: the columns of the data set which belong to the control group (in the unpaired case) or to the "before
#    treatment"-measurements
# paired: paired or unpaired data. If paired=TRUE, i.e. there are paired observations, x and y must have the same
#         length and (x[i],y[i]) will be interpreted as an observation pair
# na.rm: if na.rm=FALSE, each missing value in the data will be replaced by the genewise mean. If FALSE, the W-value
#        of every gene with one or more missing values is set to NA
# zero.rand: if there are Zeros in the paired case, what should be done with these Zeros? If zero.rand=TRUE, a sign
#            is randomly assigned to each observation pair (x[i],y[i]) with x[i]-y[i]=0. If FALSE, the sign of such
#            observation pairs is set to '-' (suggested by Lehmann (1975))
# rand: is only necessary in the paired case for the random choice of the signs of Zeros (see zero.rand)
# use.weights: if TRUE, weights will be used in the computation of p0 (for details see documentation of the
#              estimation of p0)
# delta: a vector of values for which the number of significant genes and the FDR is computed. By default,
#        these statistics are computed for any integer delta which appears in the analysis
# graphic: if TRUE, the SAM plot and the plots of of delta vs. FDR and delta vs. #significant genes is generated
# pty.square: if TRUE, a square SAM plot will be generated with x- and y-axes having the same range
# thres: the deltas for which two lines parallel to the 45°-line is generated
# use.numbers: if TRUE, the symbols for the points will be the numbers of observations which correspond to this
#              points. Otherwise circles are used. This option should only be used, if wilc=TRUE. 
# helplines: if TRUE, helplines will be generated in both the delta vs. FDR and the delta vs. #significant
#            genes plot


sam.wilc<-function(data,x,y,paired=FALSE,na.rm=FALSE,zero.rand=TRUE,rand=NA,use.weights=TRUE,delta=1:max(W.diff),
                graphic=TRUE,pty.square=TRUE,thres=round(quantile(2:max(W.diff),(0:3)/3)),use.numbers=TRUE,helplines=TRUE){
    wilc.out<-wilc.cal(data,x,y,paired=paired,zero.rand=zero.rand,rand=rand,na.rm=na.rm)
    W<-wilc.out$W        # observed unsorted W-values
    W.sort<-sort(W)      # observed sorted W-values
    W.exp<-wilc.out$W.exp   # expected (sorted) W-values
    var.0.genes<-wilc.out$var.0.genes  # index of genes with variance Zero
    n.genes<-sum(!is.na(W))    # number of genes with non missing W-value
    W.exp.number<-n.genes*wilc.out$f.null    # vector of expected numbers of the W-values under the null
    W.exp.value<-as.numeric(names(W.exp.number))    # vector of the possible W-values
    n.exp<-length(W.exp.number)  
    vec.lambda<-NULL
    vec.p0<-NULL
    for(i in 1:floor(n.exp/2)){   # estimation of p0
        vec.p0[i]<-sum(W>=W.exp.value[i] & W<=W.exp.value[n.exp-i+1],na.rm=TRUE)/sum(W.exp.number[i:(n.exp-i+1)])
        vec.lambda[i]<-1-sum(W.exp.number[i:(n.exp-i+1)])/n.genes
    }
    library(modreg)
    weights<-if(!use.weights) rep(1,length(vec.lambda)) else 1-vec.lambda
    spline.out<-smooth.spline(vec.lambda,vec.p0,w=weights,df=3)
    p0<-min(predict.smooth.spline(spline.out,1)$y,1)    
    W.diff<-W.sort-W.exp   # computation of W.diff  (we are looking for |W.diff|>=delta)
    table.W.exp<-cbind(W.exp.value,W.exp.number)
    sam.fdr.out<-sam.fdr(W.sort,W.exp,table.W.exp,p0,delta=delta,wilc=TRUE,graphic=FALSE)
    FDR<-sam.fdr.out$mat.fdr
    print(sam.fdr.out$tab.fdr)
    
    table.count<-table(W.exp,W.sort)  # for the SAM Plot the numbers of observations which correspond to the
    W.exp.value<-as.numeric(dimnames(table.count)[[1]])  # possible points is computed
    W.value<-as.numeric(dimnames(table.count)[[2]])
    n.exp<-length(W.exp.value)
    mat.count<-matrix(c(rep(W.exp.value,length(W.value)),rep(W.value,each=n.exp), 
        as.vector(table.count)),ncol = 3)
    mat.count<-mat.count[-which(mat.count[,3]==0),]
    n.sig<-NULL    
    # SAM Plot
    if(graphic){
        sam.plotter(mat.count[,2],mat.count[,1],thres,pty.square=pty.square,main="SAM Plot using Wilcoxon Rank Statistics",
            color=2:(length(thres)+1),make.legend=TRUE,wilc=TRUE,use.numbers=use.numbers,count=mat.count[,3])
        X11() 
        # delta vs. FDR and delta vs. #significant genes
        roller.coaster(FDR,helplines=helplines) # Delta vs. FDR and Delta vs. #significant genes
    }
    invisible(return(W,W.exp,W.sort,var.0.genes,W.exp.number,p0,spline.out,FDR,weights,mat.count,n.sig,
        x,y,paired,use.numbers,rand,table.W.exp,vec.p0,vec.lambda))
}
