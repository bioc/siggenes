# Copyright (C) 2002 Holger schwender

# After finding the optimal a0, this function can be used to do the Empirical Bayes Analysis of Microarray
# Experiments as proposed in Efron et al.(2001). This analyis is done for the original, untransformed data.

# a0.out: The output of find.a0
# data: the used data set; every column of this data set must correspond to one gene
# a0: if NA, the in find.a0() suggested a0 is used. Another a0 can be specified.
# p0: prior; probability that a gene is unaffected. If NA, a simple estimate of p0 will be used 
#     (min(f(Z)/f0(Z))). A better estimate for p0 can be found by using the method proposed in 
#     Efron et al.(2001) in Remark F.
# delta: if NA, the same delta will be used as in the previous analysis with find.a0(). A observation Z
#        will be called significant, if p1(Z) > delta.
# stable: if TRUE, p0 is computed by the algorithm of Storey and Tibshirani (2003a). If FALSE, the simple estimator
#         of Efron et al. (2001b) is used.
# number.int: the number of equally spaced intervals which is used in the logistic regression for the
#             calculation of the ratio f0/f. The intervals are equally spaced between min(Z) and max(Z).
# local.bin: to estimate the local FDR for Z, the proportion of the Z which fall in [Z-local.bin, Z+local.bin]
#            and the proportion of z which fall in the same interval are calculated.
# col.accession: if col.accession is a positive integer, this column of data is interpreted as the accession
#                number of the gene and it is added to the output. To avoid this, set col.accession=NA
# col.gene.name: if col.gene.name is a positive integer, this column of data is interpreted as the name of
#                the gene and is added to the output. To avoid this, set col.gene.name=NA
# q.values: if TRUE, for each gene its q-value is computed                                                
# R.fold: if TRUE, the fold change of each significant gene is calculated and added to the output
# R.dataset: the data set that is used in the computation of the R.fold. By default the data set used in the
#            other computations
# na.rm: if na.rm=FALSE, the R.fold of genes with one or more missing values will be set on NA. 
# file.out: results are stored in this file. If NA, no results will be stored in a file. 


ebam<-function(a0.out,data,a0=NA,p0=NA,delta=NA,stable=TRUE,number.int=139,local.bin=.1,col.accession=NA,
		col.gene.name=NA,q.values=TRUE,R.fold=TRUE,R.dataset=data,na.rm=FALSE,file.out=NA){
	r<-a0.out$r          # for faster calculation some of the results from previous analysis will be used
	s<-a0.out$s
	r.perm<-a0.out$r.perm
	s.perm<-a0.out$s.perm
	paired<-a0.out$paired
	if(is.na(delta))     
		delta<-a0.out$delta
	if(is.na(a0))
		a0<-a0.out$a0
	if(a0<0)             # checks if a0 is non-negative  
		stop("a0 must be larger or equal to 0.")
	if(delta>=1 || delta<=0)  # checks if delta is a probability
		stop("delta must be between 0 and 1.")
	B<-length(r.perm)/length(r)     # calculation of the number of permutations
	Z.unsorted<-r/(s+a0)         # calculation of the observed Z-values
	Z<-sort(Z.unsorted)
	z.unsorted<-as.vector(r.perm/(s.perm+a0))
	z<-sort(z.unsorted)  # calculation of the permuted z-values
	z[which(z<min(Z))]<-min(Z)    # to avoid that the logistic regression estimates will be screwed up,
	z[which(z>max(Z))]<-max(Z)    # some adjustments are made
	ratio.out<-ratio.est(Z,z,p0=p0,stable=stable,number.int=number.int)  # calculation of the posterior of the Z-values
	mat.post<-ratio.out$mat.post
	mat.repeat<-ratio.out$mat.repeat
	p0<-ratio.out$p0
	optim.out<-ratio.out$optim.out
	sig.center<-which(mat.post[,"posterior"]>delta)    # index of the "significant" centerpoints 
	nsig<-sum(mat.post[sig.center,"success"])          # calculation of the number of significant genes and...
	false<-sum(mat.repeat[sig.center,"n"]-mat.repeat[sig.center,"success"])/B # ... of the number of falsely
	                                                                          # called genes
	fdr<-p0*false/nsig   # calculation of the FDR
	posterior<-rep(mat.post[,"posterior"],mat.post[,"success"])  # calculation of the posterior for each
	mat.post.Z<-cbind(Z,posterior)                               # significant gene
	sig.genes<-which(mat.post.Z[,"posterior"]>=delta)     # index of the significant (sorted) genes
	FDR<-c(p0=p0,nsig=nsig,false=false,fdr=fdr)   # for output
	local2<-NULL            # estimation of the local FDR by binning the observations into small intervals
	for(i in 1:length(sig.genes)){
		local2[i]<-p0*sum(z>=Z[sig.genes[i]]-local.bin & z<=Z[sig.genes[i]]+local.bin)/
			(B*sum(Z>=Z[sig.genes[i]]-local.bin & Z<=Z[sig.genes[i]]+local.bin))
	}
	# output is made
	mat.ebam<-cbind("Z"=round(mat.post,Z[sig.genes,1],4),"p1(Z)"=round(mat.post.Z[sig.genes,2],4),
		"local1"=round(1-mat.post.Z[sig.genes,2],4),"local2"=round(local2,4))
	if(!is.na(col.accession))
		mat.ebam<-cbind("access"=data[row.sig.genes,col.accessiom],mat.ebam)
	if(q.values){
		mat.qvalue<-q.value.cal(Z.unsorted,z.unsorted,p0)$mat.qvalue
		q.value<-mat.qvalue[order(mat.qvalue[,1]),3]
		mat.ebam<-cbind(mat.ebam,"q-value"=round(q.value[sig.genes],5))
	}
	if(R.fold){
		fold.change<-R.fold.cal(R.dataset[row.sig.genes,],a0.out$x,a0.out$y,na.rm=na.rm)
		mat.ebam<-cbind(mat.ebam,"R-fold"=round(fold.change[,3],4))
	}
	if(!is.na(col.gene.name))
		mat.ebam<-cbind(mat.ebam,"gene"=substring(data[row.sig.genes,col.gene.name],1,50))
	mat.ebam<-cbind("ID"=row.sig.genes,mat.ebam)
	ebam.out<-as.data.frame(mat.ebam)    
	names(ebam.out)<-c("accession","Z","posterior","local1","local2")
	cat("Using a0 =",round(a0,4),"and the original Z values, there are",nsig,"significant genes and",false,
		"falsely called genes.","\n","For p0 =",round(p0,4),", the FDR is",round(fdr,4),".","\n")
	if(!is.na(file.out)){   # output is stored in a file
		cat("Results of the Empirical Bayes Analysis of Microarray Experiments","\n","\n","\n",
			"Significance criterion: p1(Z) >=",delta,"\n","\n","a0:",round(a0,4),"\n","p0:",round(p0,4),"\n",
			"significant genes:",nsig,"\n","falsely called genes:",round(false,4),"\n","FDR:",round(fdr,4),
			"\n","\n","\n","Genes called significant:","\n","\n",file=file.out)
		write.table(ebam.out,file=file.out,sep="\t",append=TRUE,row.names=FALSE,col.names=TRUE,quote=FALSE)
		cat("\n","\n","local1: estimation of local FDR using the logistic regression estimates",
			"\n","local2: local FDR is estimated by binning the observations into small intervals",
			file=file.out,append=TRUE)
		cat("\n","\n","Output is stored in",file.out,"\n")
	}
	plot(mat.post.Z[,"Z"],mat.post.Z[,"posterior"],main="Posterior probability for Z values",
		xlab="Z values",ylab="p1(Z)")     # the posterior of the genes is plotted
	abline(h=delta,lty=2)
	# mark the significant genes with green color
	points(mat.post.Z[which(mat.post.Z[,2]>delta),1],mat.post.Z[which(mat.post.Z[,2]>=delta),2],col=4)
	mat.Z.unsorted<-cbind(Z.unsorted,mat.post.Z[rank(Z.unsorted),2])
	invisible(return(list(mat.repeat=mat.repeat,optim.out=optim.out,
                              mat.post.Z=mat.post.Z,ebam.out=ebam.out,
                              FDR=FDR,a0=a0,mat.Z.unsorted=mat.Z.unsorted,
                              Z.unsorted=Z.unsorted,row.sig.genes=row.sig.genes,p0=p0)))
}
