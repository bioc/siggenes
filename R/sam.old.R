# Copyright (c) 2002 Holger Schwender, University of Dortmund, Germany

# This function does the Significance Analysis of Microarray Experiments. The output of this function is
# a table of statistics (p0, #number of significant genes, #number of falsely called genes, #FDR) for some
# delta and a SAM Plot for some (not necessarily the same) delta. Additionally it is possible that this
# function computes a delta for a given number or proportion of genes which should be significant.

# The output of this function can be used to select the delta in the SAM Analysis. For further analysis
# the output of sam() must be assigned to an object. 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# data: the data set; there is only one condition: every row of this data set must represent a gene
# x: the columns of the data set which belong to the cases (in the unpaired case) or the "after treatment"- 
#    values (in the paired case) 
# y: the columns of the data set which belong to the control group (in the unpaired case) or to the "before
#    treatment"-measurements
# paired: paired or unpaired data
# mat.samp: the permutation matrix. If specified and correct, this matrix will be used, even if rand and B are
#           specified.
# B: number of permutations used in the calculation of the Null. Will not be used, if mat.samp is specified.
# balanced: if TRUE, balanced permutations will be used
# na.rm: if na.rm=FALSE, the d-values of genes with one or more missing values will be set to NA. If na.rm=TRUE, the
#        missing values will be removed during the computation of the d-values.
# s0: the fudge factor. If NA, s0 will automatically be computed as the minimum coefficient of variation of
#     d(i) as a function of s(i)
# *.s0: These options are necessary for the calculation of s0 in fudge(). For detailed information see fudge().
#       alpha.s0 are the possible values of s0; if include.s0=TRUE, then s0=0 is also a possible choice;
#       factor.s0 is the constant with which the MAD is multiplied
# p0: the (prior) probability that a gene is unaffected
# *.p0: These options are used in p0.est(). For detailed information see p0.est().
#       lambda.p0 is a number between 0 and 1 (including 0 and 1) which is used to estimate p0; if lambda.p0=1,
#       a natural cubic spline with 3 df of p0(vec.lambda.p0[i]) on vec.lambda.p0[i] is fitted
# *.fdr: These options are used in sam.fdr(). For detailed information see sam.fdr().
#        delta.fdr are the delta for which the FDR are computed; if med.fdr=TRUE, the median number, otherwise
#        the expected number, of falsely called genes is calculated; if graphic.fdr=TRUE, the plots will be made;
#    for each thres.fdr two lines parallel to the 45°-line are plotted; if pty.fdr=TRUE, a square SAM plot
#    is made; if help.fdr=TRUE, helplines will be plotted in both the delta vs. FDR and the delta vs.
#    #significant genes plot
#    If you only like to use ngenes and avoid the computation of the FDR for some delta, i.e the analysis
#        which is done by sam.fdr(), set delta.fdr=NULL
# ngenes: a number or proportion of genes for which a delta is computed so that about this number or percentage
#         of genes fall outside the thresholds cutup and cutlow, i.e. so that about this number or proportion
#         of genes is called significant. Default is NA, i.e. no such analysis will be done.
#           This is a option in sam.ngenes().  
# iteration: the number of iterations used in sam.ngenes(). For details see sam.ngenes().
# initial.delta: the initial guesses for delta in sam.ngenes(). For details see sam.ngenes()
# rand: the set.seed. Default is NA, i.e. there will be no set.seed().      


sam.old<-function(data,x,y=NULL,paired=FALSE,mat.samp=NULL,B=100,balanced=FALSE,na.rm=FALSE,s0=NA,alpha.s0=seq(0,1,.05),include.s0=TRUE,
        factor.s0=1.4826,p0=NA,lambda.p0=1,vec.lambda.p0=(0:95)/100, delta.fdr=(1:10)/5,
        med.fdr=TRUE,graphic.fdr=TRUE,thres.fdr=seq(0.5,2,.5),pty.fdr=TRUE,help.fdr=TRUE,ngenes=NA,iteration=3,
        initial.delta=c(0.1,seq(.2,2,.2),4),rand=NA){
    rs.out<-rs.cal(data,x,y,paired=paired,mat.samp=mat.samp,B=B,bal=balanced,na.rm=na.rm,rand=rand)
    r<-rs.out$r     # calculation of the observed and expected r- and s-values
    s<-rs.out$s
    r.perm<-rs.out$r.perm
    s.perm<-rs.out$s.perm
    mat.samp<-rs.out$mat.samp   # the permutation matrix
    var.0.genes<-rs.out$var.0.genes  # the row numbers corresponding to the genes with variance Zero
    if(is.na(s0))   # calculation of the fudge factor s0
        s0<-fudge(r,s,alpha=alpha.s0,include.zero=include.s0,factor=factor.s0)$s.zero
    d<-r/(s+s0)      # observed unsorted d-values
    if(any(is.na(d)))
        cat("There are",sum(is.na(d)),"missing d values.","\n","\n")
    d.sort<-sort(d,na.last=TRUE)  # sorted observed d-values
    d.perm<-r.perm/(s.perm+s0)  
    d.perm<-apply(d.perm,2,sort,na.last=TRUE)  # matrix of sorted permuted d-values
    d.bar<-rowMeans(d.perm)   # expected sorted d-values
    if(is.na(p0))
        p0<-p0.est(d,d.perm,lambda=lambda.p0,vec.lambda=vec.lambda.p0)$p0  # estimation of p0
    
    FDR<-NULL
    delta.ngenes<-NULL
    fdr.ngenes<-NULL
    if(!is.null(delta.fdr)){   # doing a SAM Analysis for some delta
        cat("SAM Analysis for a set of delta:","\n")
        sam.fdr.out<-sam.fdr(d.sort,d.bar,d.perm,p0,delta=delta.fdr,med=med.fdr,graphic=graphic.fdr,thres=thres.fdr,
            pty.square=pty.fdr,helplines=help.fdr)
        FDR<-sam.fdr.out$mat.fdr    # calculation of some statistics like #significant, #falsely called 
                            # genes, FDR
        print(sam.fdr.out$tab.fdr)  #output
        cat("\n")
    }
    if(!is.na(ngenes)){   # for given number or percentage of significant genes a delta is computed for the
                              # calculation of the FDR
        sam.ngenes.out<-sam.ngenes(d.sort,d.bar,d.perm,ngenes=ngenes,iteration=iteration,initial.delta=initial.delta,med=med.fdr,
            p0=p0)
            fdr.ngenes<-sam.ngenes.out$fdr.ngenes       # some statistics for the exact delta or for the upper
          			                        # and lower bound
            delta.ngenes<-sam.ngenes.out$delta.ngenes   # exact delta
          
    }
    if(is.null(FDR))
        FDR<-fdr.ngenes
    structure(list(d=d,d.sort=d.sort,s=s,d.bar=d.bar,d.perm=d.perm,mat.samp=mat.samp,s0=s0,
		FDR=FDR,p0=p0,fdr.ngenes=fdr.ngenes,delta.ngenes=delta.ngenes,med.fdr=med.fdr,
		x=x,y=y,paired=paired,var.0.genes=var.0.genes))
}
        
