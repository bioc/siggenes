# Copyright (C) 2002 Holger Schwender

# To get the same results in sam() in S and R, one has to use the same permutation matrix.
# Unfortunately the random numbers generators in S and R work differently.
# This program solves this problem by defining the permutation matrix 'mat.samp' which can
# be used in sam.bal() to do a SAM analysis. The difference between sam() and sam.bal() is
# that in sam.bal() a predefined permutation matrix can be used where in sam() this matrix
# will be calculated.

# CAUTION: The results in S-plus differ from the results in R - even for the same set.seed().

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# n.x: the number of cases (i.e. length of 'case' in sam())
# n.y: the number of controls ( i.e. length of 'control' in sam())
# B: the number of permutations which should be used in the SAM analysis
# paired: if paired data or not
# balanced: if balanced is TRUE, balanced permutations or permutations which are as balanced as possible
#           are used
# rand: here one can define set.seed(), if NA, set.seed() will not be used
# file.out: the permutation matrix will be stored here for further use e.g. in R 
#           if there shouldn't be any output, set file.out=NA 

sam.sampler<-function(n.x,n.y,B,paired=FALSE,balanced=FALSE,rand=NA,file.out=NA){
    if(!is.na(rand))
        set.seed(rand)
    if(!paired){
        if(!balanced){
            mat.samp<-matrix(0,B,n.x+n.y)
            for(i in 1:B)
                mat.samp[i,]<-sample(c(rep(1,n.x),rep(0,n.y)),n.x+n.y)
        }
        if(balanced){
            mat.samp.x<-matrix(0,B,n.x)
            mat.samp.y<-matrix(0,B,n.y)
            for(i in 1:B){
                mat.samp.x[i,]<-sample(rep(c(1,0),ceiling(n.x/2)),n.x)
                mat.samp.y[i,]<-sample(rep(c(1,0),ceiling(n.y/2)),n.y)
            }
            mat.samp<-cbind(mat.samp.x,mat.samp.y)
        }
    }
    if(paired){
        if(n.x!=n.y)
            stop("x must be equal to y")
        mat.samp<-matrix(0,B,n.x)
        for(i in 1:B)
            mat.samp[i,]<-sample(c(-1,1),n.x,replace=TRUE)
    }
    if(!is.na(file.out))
        write.table(mat.samp,file=file.out,sep="\t")
    return(mat.samp)
}
    
