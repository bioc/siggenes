# Copyright (c) 2002 Holger Schwender

# This function estimates p0, the probability of an unaffected gene, as described in Storey (2002)
# "False Discovery Rates and Q-values for Inference in DNA Microarray Experiments"

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# d: the observed d-values
# d.perm: matrix of the permuted d-values
# lambda: number between 0 and 1. p0 is computed as a function of lambda. If lambda=1, p0(1) is computed
#         using natural cubic splines. This computation leads to an optimal estimation of p0.
# vec.lambda: if lambda=1, a natural cubic spline with 3 df of p0(vec.lambda[i]) on vec.lambda[i] is fitted 


p0.est<-function(d,d.perm,lambda=1,vec.lambda=(0:95)/100){
	if(lambda>1 || lambda<0)   # limitation of lambda
		stop("lambda has to be between 0 and 1.")
	if(lambda!=1)       # if lambda!=1, one is only interested in p0(lambda)
		vec.lambda<-lambda
	vec.p0<-NULL
	d<-na.exclude(d)
	m<-length(d)
	quan<-quantiles(na.exclude(d.perm),c(vec.lambda/2,1-rev(vec.lambda)/2))  # calculation of the null quantiles
	for(i in 1:length(vec.lambda))
		vec.p0[i]<-sum(d>quan[i] & d<quan[length(quan)-i+1])/((1-vec.lambda[i])*m)   # calculation of p0(lambda)
	if(lambda!=1){
		p0<-min(vec.p0,1)      # if lambda!=1, one is only interested in p0(lambda); truncate p0 at 1
		return(list(p0=p0,vec.p0=vec.p0))
	}
	spline.out<-smooth.spline(vec.lambda,vec.p0,w=1-vec.lambda,df=3)   # smooth natural cubic splines with 3 df
	p0<-min(predict(spline.out,1)$y,1)   # compute p0(1) which is the estimation of p0
	structure(list(p0=p0,spline.out=spline.out,vec.p0=vec.p0))
}
