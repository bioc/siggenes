# Copyright (C) 2002 Holger schwender

# This program calculates the logistic regression estimates of the ratio f0/f as suggested
# by Efron et al. (2001).

# This program is be used by ebam().

# Required libraries/functions: library(splines), logLik.repeat

# Z.norm: the vector of the normal score transformed observed Z-values, i.e. f(Z.norm) is a nearly perfect
#     N(0,1) density.
# z.norm: the vector of the transformed null scores z, using the same transformation as above.
# p0: prior; probability that a gene is unaffected. If not specified the usual estimation (see below) will
#     be done.
# stable: if TRUE, p0 is estimated by the algorithm of Storey an Tibshirani (2003a). If FALSE, the simple estimate
#         of Efron et al. (2001b) is used.
# number.int: The number of intervals intervals which should be used to aggregate the data. Default is 139,
#             the number of intervals Efron et al. (2001) used.
       


ratio.est<-function(Z.norm,z.norm,p0=NA,stable=TRUE,number.int=139){
    library(splines)   # needed for ns()
    min.int<-floor(100*min(Z.norm))/100     # lower limit for the first interval
    max.int<-ceiling(100*max(Z.norm))/100   # upper bound for the last interval
    #if(method=="e")      # method of equal intervals
        interval<-seq(min.int,max.int,length=number.int+1)      # calculation of the intervals
    #if(method=="n"){     # method of equal number of genes per interval
    #   index<-intervaller(length(Z.norm),number.int)    # genes whose values are the limits of the intervals
    #   interval<-c(min.int,Z.norm[index])      # calculation of the interval borders
    #   interval[length(interval)]<-max.int     # set the upper bound for the last interval to max.int
    #}
    # calculation of the centerpoints
    center<-(interval[2]-interval[1])/2+interval[-length(interval)]
    bin.Z<-cut(Z.norm,interval,include.lowest=TRUE)       # bin the Z values into the intervals
    bin.z<-cut(z.norm,interval,include.lowest=TRUE)       # bin the z values into the intervals
    success<-tabulate(bin.Z,length(levels(bin.Z)))     # for each interval get the number of Z and z values
    failure<-tabulate(bin.z,length(levels(bin.z)))     # which belong to this interval. Z is a success, z is
                                                       # a failure
    n<-success+failure     # number of values for each interval
    p<-success/n           # proportion of Z values for each interval
    log.bino<-lgamma(n+1)-(lgamma(success+1)+lgamma(n-success+1)) # calculation of the logarithm of the binomial
                                                # coefficient in the loglikelihood; no further use
    ns.out<-ns(center,5)        # get the natural splines matrix for the centerpoints
    
    # make a data frame which will be used by ms(); exclude rows with NA, i.e. rows which belong to intervals
    # with neither a Z nor a z value
    mat.repeat<-na.exclude(cbind(log.bino,ns.out,n,success,p,center))  
    mat.repeat<-as.data.frame(mat.repeat)                              
    names(mat.repeat)<-c("log.bino","x1","x2","x3","x4","x5","n","success","p","center")
    
    # minimize the negative loglikelihood; set the start values of the parameters to 0
    attach(mat.repeat)        # method="BFGS" seems to do the best job and seems to lead to
                  # almost the same results as in S-plus 
    optim.out<-optim(rep(0,6),neglogLik.repeat,method="BFGS")
    b<-as.vector(optim.out$par)        # get the estimated parameters
    mat.model<-as.matrix(cbind(1,mat.repeat[,2:6]))  # get the model matrix
    pi.Z<-exp(mat.model%*%b)/(1+exp(mat.model%*%b))    # calculate pi(Z); the probability of a success in an
                                                       # interval
    B<-length(z.norm)/length(Z.norm)     # calculation of the number of permutation used to estimate the Null
    
    if(is.na(p0)){
        if(stable)
            p0<-p0.est(Z.norm,z.norm)$p0
        else
            p0<-min((B*pi.Z)/(1-pi.Z))
        }
    posterior<-1-p0*(1-pi.Z)/(B*pi.Z)     # calculate the posterior probability 
    posterior[which(posterior<0)]<-0      # truncate posterior at 0
    mat.post<-cbind(mat.repeat[,c("center","success")],posterior)
    structure(list(mat.repeat=mat.repeat,optim.out=optim.out,p0=p0,mat.post=mat.post))
}
