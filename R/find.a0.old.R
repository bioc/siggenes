# Copyright (C) 2002 Holger Schwender

# This function is looking for the optimal fudge factor a0. Efron et al.(2001) defines that the optimal choice
# is the a0 which leads to the most significant genes, i.e. most genes with a posterior probability p1(Z)>0.9.
# a0 is either set to 0 or to a quantile of the s-values, i.e. the standard deviations of the genes.

# find.a0 suggest a choice of a0 for further analysis based on the above optimalization criterion. But one
# can choose another value of a0 in further analysis. For confirmation the logit of the posterior probabilities
# is plotted.

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
# delta: each observation with posterior probability p1(Z) >= delta is called significant. Default is 0.9 which
#        is used in Efron et al.(2001)
# alpha: these alpha quantiles of the s-values will be used to find the optimal a0
# include.0: if TRUE, a0 = 0 is also considered in the search for an optimal a0
# p0: prior; the probability that a gene is unaffected. If not specified, only a very simple estimation for
#     p0 is used, i.e. p0 = min(f(z)/f0(z))  
# stable: if TRUE, p0 is computed by the algorithm of Storey and Tibshirani (2003a). If FALSE, the simple estimator
#         of Efron et al. (2001b) is used.
# number.int: number of equally spaced intervals (between min(Z) and max(Z), where Z are the observed values
#             of the genes) which are used in the logistic regression estimate of the ratio f0/f.
#             Default is 139, which is used in Efron et al.(2001)
# rand: specifies the set.seed for the calculation of the permutation matrix. Will not be used, if mat.samp
#       is specified
# plot.legend: if TRUE, there will be a legend in the logit(posterior)-plot. It is highly recommended to set
#              plot.legend to FALSE, if length(alpha) exceeds 10. 




find.a0.old<-function(data,x,y,paired=FALSE,mat.samp=NULL,B=100,balanced=FALSE,na.rm=FALSE,delta=0.9,alpha=(0:9)/10,include.0=TRUE,p0=NA,stable=TRUE,number.int=139,rand=NA,plot.legend=TRUE){
    rs.cal(data,x,y,paired=paired,mat.samp=mat.samp,bal=balanced,B=B,na.rm=na.rm,rand=rand)->rs.out
    r<-rs.out$r        # The r- and s-values for the numerator and the denominator, respectively, are calculated
    s<-rs.out$s        # for both the observed and the permuted observations.
    r.perm<-rs.out$r.perm
    s.perm<-rs.out$s.perm
    mat.samp<-rs.out$mat.samp
    vec.a0<-quantile(s,alpha,na.rm=TRUE) # The vector with the values of a0, which are used to find the optimal a0, is built
    if(include.0)
        vec.a0<-c(0,vec.a0)
    n.genes<-length(na.exclude(r))      # number of genes with non-missing values
    sig.a0<-NULL
    mat.post<-NULL
    Z.norm<-qnorm(((1:n.genes)-3/8)/(n.genes+0.25))  # normal score transformation using Blom normal scores
    for(i in 1:length(vec.a0)){     # calculation of the posterior for the different a0
        Z<-sort(r/(s+vec.a0[i]))         # calculation of the observed Z-values
        z<-sort(as.vector(r.perm/(s.perm+vec.a0[i])))  # calculation of the permuted z-values
        z.norm<-approx(Z,Z.norm,z,rule=2)$y   # use linear interpolation to do the transformation for the
                            # z-values which was used to transform the Z-values
                            # rule=2 means that every z which is smaller than min(Z) is set to min(Z.norm)
                            # and every z which is larger than max(Z) is set to max(Z.norm)
        mat.ratio<-ratio.est(Z.norm,z.norm,p0=p0,number.int=number.int)$mat.post  # posterior is calculated
        mat.post<-rbind(mat.post,cbind(a0=rep(vec.a0[i],nrow(mat.ratio)),mat.ratio))  
        sig.a0[i]<-sum(mat.ratio[which(mat.ratio[,"posterior"]>=delta),"success"]) # number of significant genes
                                                   # using the current a0 is calculated
    }
    logit.post<-log(mat.post[,"posterior"]/(1-mat.post[,"posterior"]))   
                # logit(posterior) will be used in the plot to emphasize the differences in the tails
    plot(mat.post[which(mat.post[,"a0"]==vec.a0[1]),"center"],logit.post[which(mat.post[,"a0"]==vec.a0[1])],
        main="Transformed Z values vs. Logit of the Posterior",xlab="transformed Z values",
        ylab="logit(posterior)",type="l",xlim=c(-4,4),ylim=c(0,max(logit.post[which(logit.post!=Inf)])+0.5))
    if(any(logit.post==Inf))
	cat("Warning: Some of the logit posterior probabilities are Inf. These probabilities are not plotted.",
		"\n","\n")
    for(i in 2:length(vec.a0))       # logit posterior is plotted
        lines(mat.post[which(mat.post[,"a0"]==vec.a0[i]),"center"],logit.post[which(mat.post[,"a0"]==vec.a0[i])],
            col=i)
    abline(h=log(delta/(1-delta)),lty=4)    # this line corresponds to 0.9 in p1(Z) >= 0.9
    vec.a0.names<-NULL
    vec.a0.name<-NULL
    for(i in 1:length(alpha)){     # "nice" names for the a0 are made
        vec.a0.name[i]<-paste(c("a0=",round(vec.a0[i+1],4)," (alpha=",alpha[i],")"),collapse="")
        vec.a0.names[i]<-paste(c("alpha=",alpha[i]," (",sig.a0[ifelse(include.0,i+1,i)],")"),collapse="")
    }
    if(include.0){
        vec.a0.name<-c("a0=0",vec.a0.name)
        vec.a0.names<-c(paste(c("a0=0 (",sig.a0[1],")"),collapse=""),vec.a0.names)
    }
    if(plot.legend)      # a legend is plotted; highly recommended: set to FALSE, if length(alpha)>10
        legend(-1.1,max(logit.post)+.6,legend=vec.a0.names,lty=1,cex=0.8,col=1:length(vec.a0),bty="n")
    names(sig.a0)<-vec.a0.name
    cat("\n","Number of significant genes for some a0:","\n")   # the most important information is displayed
    print(sig.a0)
    a0<-vec.a0[which(sig.a0==max(sig.a0))][1]
    cat("\n","Suggested choice for a0:",round(a0,4))         # a choice for a0 is suggested
    if(a0!=0)
        cat("   (the",names(a0),"quantile of the s-values)")  
    cat("\n")
    structure(list(r=r,s=s,r.perm=r.perm,s.perm=s.perm,mat.samp=mat.samp,sig.a0=sig.a0,a0=a0,
		delta=delta,vec.a0=vec.a0,x=x,y=y))
}
        
