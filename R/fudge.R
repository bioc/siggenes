# Copyright (c) 2002 Holger Schwender

# This program calculates the fudge factor s0 for further use in calculation of the d-values. Adding s0 to
# rhe denominator of the d-values ensures that the variance of d(i) is independent of gene expressions.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# r: a vector of r-values. For example, for the two class unpaired data r(i) = mean(expression value of
#    gene i in class 1) - mean(value of gene i in class 2). r(i) is the numerator of d(i)
# s: a vector of s-values. s(i) is the standard deviation of gene i. s(i)+s0 is the denominator of d(i)
# alpha: these alpha quantiles will be used (as values for s0) to find the optimal s0
# include.zero: if TRUE, s0=0 is also considered in the search for an optimal a0
# factor: number with which the median of the absolute values is multiplied. The MAD is used in the 
#         calculation of the coefficient of variance of d(i) as a function of s(i) in this analysis.
#         The default value makes the estimate consistent for the standard deviation at the Gaussian model.


fudge<-function(r,s,alpha=seq(0,1,0.05),include.zero=TRUE,factor=1.4826){
if(max(alpha)>1 || min(alpha)<0){  
    stop("alpha has to be between 0 and 1") 
}
if(any(round(100*alpha,10)!=round(100*alpha,0))){    # alpha has to be a percentile   
    cat("Warning: At least one alpha is not a percentile. Only the first two decimal digits are retained.","\n")
            alpha<-signif(alpha,2)
}
quan<-quantile(s,seq(0,1,0.01),na.rm=TRUE)   # the percentiles of s are calculated
cv<-NULL
for(i in 1:length(alpha)){  # for the alpha quantile of the s-values the coefficient of variation is calculated
    v<-NULL
    for(j in 1:100){       # compute the MAD of the d(i) as a function of the s(i)
        d.alpha<-r[which(s>=quan[j] & s<quan[j+1])]/(s[which(s>=quan[j] & s<quan[j+1])]+quan[100*alpha[i]+1])
        v[j]<-mad(d.alpha,constant=factor) 
    }   
    cv[i]<-sqrt(var(v))/mean(v)  # compute the coefficient of variation of these MAD values
}
if(include.zero){  # the same for s0=0
    v<-NULL
    for(j in 1:100){
        d.alpha<-r[which(s>=quan[j] & s<quan[j+1])]/s[which(s>=quan[j] & s<quan[j+1])]
        v[j]<-mad(d.alpha,constant=factor)}
    cv.zero<-sqrt(var(v))/mean(v)
    if(cv.zero<min(cv)){  # is s0=0 the best choice?
        cat("s0 =",0,"\n","\n")  # some output things
        s.zero<-0
        alpha.hat<-NA
}}
if(!include.zero || cv.zero>=min(cv)){  # again some output
    alpha.hat<-alpha[which(cv==min(cv))]     # which alpha quantile of the s-values is the best choice for s0?
    s.zero<-quan[100*alpha.hat+1]
    cat("alpha.hat: ",alpha.hat,"\n","      s0: ",round(s.zero,4),"\n","\n")
}
if(!include.zero)
    cv.zero<-NA
invisible(return(alpha.hat,s.zero,cv,cv.zero))
}
