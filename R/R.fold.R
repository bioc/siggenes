# Copyright (c) 2002 Holger Schwender

# This function computes the fold change. Here the R-fold criterion is only used as an additional criterion.
# The fold change is reported but a gene will not be called not-significant if this gene does not fullfil
# the criterion that |mean(x_i1)/mean(x_i2)|>=t or <=t for some t, where x_i1 and x_i2, respectively, are
# the average expression levels of gene i under two different conditions.

# The fold change for a gene i will be set to NA, if either mean(x_i1) or mean(x_i2) is less or equal to 0
# because such fold changes cannot be unambiguously computed.

# data: the used data set; this data set must be the same data set as in sam(), but it could be, e.g.,
#        the unnormalized version of the data set used in sam() if this data set was normalized
# x: the vector of columns which belong to the cases (unpaired) or "after treatment"-measurements (paired)
# y: the vector which contains the columns that belong to the controls (unpaired) or to the "before
#    treatment"-measurements (paired)
# na.rm: if na.rm=FALSE, the d-values of genes with one or more missing values will be set to NA. If na.rm=TRUE, the
#        missing values will be removed during the computation of the d-values.


R.fold.cal<-function(data,x,y,na.rm=FALSE){
    X<-as.matrix(data[,c(x,y)])  # for easier calculation
    mode(X)<-"numeric"
    mean.x <- rowMeans(X[,1:length(x)],na.rm=na.rm)  # compute the genewise mean of the first group
    mean.y <- rowMeans(X[,(length(x)+1):ncol(X)], na.rm=na.rm) # compute the genewise mean of the second group
    vec.R.fold<-mean.x/mean.y  # compute the fold changes
    vec.R.fold[which(mean.x<=0 | mean.y<=0)]<-NA   # set the fold changes to 0 which have at least
                                                            # one non-positive group mean
    mat.R.fold<-cbind(mean.x=mean.x,mean.y=mean.y,R.fold=vec.R.fold)  # for output
    invisible(return(mat.R.fold))
}
