# Copyright (c) 2002 Holger Schwender

# This function computes the q-values which were introduced by John Storey. Q-values can be 
# interpreted as a p-value in a multiple testing analysis.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# d: vector of the observed d-values
# d.perm: matrix of the permuted d-values
# p0: the probability that a gene is unaffected

q.value.cal<-function(d,d.perm,p0){
	n.with.NA<-length(d)  # number of genes
	d<-na.exclude(d)            # remove the NAs 
	d.perm<-na.exclude(d.perm)
	m<-length(d)        # number of genes
	B<-length(d.perm)/length(d)     # number of permutations
	mat.addup<-matrix(c(abs(d),abs(d.perm),rep(1,m),rep(0,m*B),sign(d),sign(d.perm)),m*(B+1),3)
	mat.addup<-mat.addup[order(mat.addup[,1]),1:3]  # a matrix is build with the absolut sorted
	                                                # observed and permuted d-values
	vec.false<-m-(which(mat.addup[,2]==1)-(1:m))/B  # the number of falsely called genes is computed
	                                       # instead of the p-values which were suggested by Storey (Nov 2002)
	q.value<-p0*vec.false[1]/m  # calculation of the q-value
	for(i in 2:m)
		q.value[i]<-min(p0*vec.false[i]/(m-i+1),q.value[i-1])
	mat.qvalue<-cbind(d=d[order(abs(d))],false=vec.false,q.value=q.value) # for further calculation and output
	if(n.with.NA!=m)  # every gene with a missing d-value gets a q-value which is NA
		mat.qvalue<-rbind(mat.qvalue,matrix(NA,n.with.NA-m,3))
	return(mat.addup,vec.false,q.value,mat.qvalue)
}
	