# Copyright (c) 2003 Holger Schwender

# If the output of a previous analysis with SAM was assigned to an object, this function can
# be used to compute the number of significant genes and the FDR for another set of thresholds
# delta.


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# sam.out: the object to which a previous analysis with sam() was assigned
# delta: a vector containing values for the threshold delta

sam.delta<-function(sam.out,delta){
	sam.fdr(sam.out$d.sort,sam.out$d.bar,sam.out$d.perm,sam.out$p0,delta=delta,med=sam.out$med.fdr,graphic=FALSE)$tab.fdr
}