# Copyright (c) 2002 Holger Schwender

# This program is required for maximize the likelihood of the logistic regression with 
# repeated observations as described in Neter et al. (1996).

# Caution: 1. This function contains the negative likelihood. 
#          2. This function only contains the for the maximization of the likelihood 
#	      important terms of the likelihood, i.e. the log(binomial coefficient) was omitted 

# b: vector of parameters for which we like to get ML-estimations


neglogLik.repeat<-function(b){sum(-success*(b[1]+b[2]*x1+b[3]*x2+b[4]*x3+b[5]*x4+b[6]*x5)
	+n*log(1+exp(b[1]+b[2]*x1+b[3]*x2+b[4]*x3+b[5]*x4+b[6]*x5)))}