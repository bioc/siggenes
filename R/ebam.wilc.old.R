# Copyright (C) 2002 Holger Schwender

# This program does an Empirical Bayes Analysis of Microarray Experiments using Wilcoxon Scores as described
# in Efron et al.(2001a), "Microarrays, Empirical Bayes Methods, and False Discovery Rate"

# data: the used data set; condition: every row of the data set must represent a gene
# x: the columns which belong to the cases (unpaired) or to the "after treatment"-measurement (unpaired)
# y: the columns which belong to the controls (unpaired) or to the "before treatment"-measurements (paired)
#    In the paired case x and y must have the same length; (x[i], y[i]) belong to each other
# delta: a gene Z will be called significant if the posterior probability p1(Z) >= delta. Default is 0.9, the
#        value Efron et al. used in their analysis
# p0: prior; probability that a gene is unaffected. If NA, a simple or a more stable estimation of p0 will be
#     calculated (see stable.p0)
# stable.p0: If TRUE, the algorithm of Storey and Tishirani (2003a) is used. If FALSE, the estimate of pi0 is computed
#            that ensures that the posterior probability of being differenntially expressed is always positive.
# use.offset: if TRUE, a Poisson regression with offset will be done. If FALSE, a Poisson regression without offset
#             is performed.
# use.weights: if TRUE, weights will be used in the computation of p0 (for details see documentation of the
#              estimation of p0)
# ties.rand: A problem of ranking is that there could be ties. If there are ties, the Wilcon test statistic W
#            could be a non-integer. If ties.rand=T, then such a statistic will be randomly assigned to either
#            floor(W) or ceiling(W). If F, the statistic is treated in a conservative way, i.e. it will be
#            assigned to ceiling(W) if W < mean of the Null and to floor(W) if W > mean of the Null.
#            To use non-integer Ws in further analysis would totally screw up the density estimation of the
#            observed Wilcoxon test statistics.
# zero.rand: Using paired data it could happen that x[i]=y[i]. If zero.rand=T, the sign of such a pair will be
#            randomly assigned for the Wilcoxon Sign-Rank Test. If F, the method of Lehmann (1975) is used.
# ns.df: the number of df for the natural splines. If ns.df is n, n-1 knots will be used in the calculation
#        of the natural splines. Default is 5, which is used by Efron et al.(2001a).
# col.accession: if col.accession is a positive integer, this column of data is interpreted as the accession
#                number of the gene and it is added to the output. To avoid this, set col.accession=NA
# col.gene.name: if col.gene.name is a positive integer, this column of data is interpreted as the name of
#                the gene and is added to the output. To avoid this, set col.gene.name=NA
# R.fold: if TRUE, the fold change of each significant gene is calculated and added to the output
# R.dataset: the data set that is used in the computation of the R.fold. By default the data set used in the
#            other computations
# file.out: some results are stored in this file. If NA, no such storage will happen
# rand: the set.seed. If NA, no set.seed will be used
# na.rm: if na.rm=FALSE, the R-fold of genes with one or more missing values will be set to NA. If na.rm=T, the
#        missing values will be removed during the computation of the R-fold.


ebam.wilc.old<-function(data,x,y,paired=FALSE,delta=.9,p0=NA,stable.p0=TRUE,use.offset=TRUE,use.weights=TRUE,ties.rand=TRUE,
		zero.rand=TRUE,ns.df=5,col.accession=NA,col.gene.name=NA,R.fold=TRUE,R.dataset=data,file.out=NA,rand=NA,na.rm=FALSE){
	library(splines)  # necessary for ns()
	if(!is.na(rand))
		set.seed(rand)
	Y<-as.matrix(data[,c(x,y)])   # for easier calculations a matrix of the data is made
	mode(Y)<-"numeric"
	n.genes<-nrow(Y)          # number of genes
	if(!paired){    # unpaired case
		n.x<-length(x)
		n.y<-length(y)
		W.mean<-n.x*(n.x+n.y+1)/2    # some statistics of the null density are calculated
		W.min<-n.x*(n.x+1)/2
		W.max<-n.x*(2*n.y+n.x+1)/2
		Y.rank<-t(apply(Y,1,rank))         # rank the observations for each gene
		y.wilk <- rowSums(Y.rank[,1:n.x])  # calculation of the Wilcoxon test statistic
		f.null<-dwilcox(0:(n.x*n.y),n.x,n.y) # calculation of the null density
	}
	if(paired){  # paired case
		n<-length(x)
		Y<-Y[,1:n]-Y[,(n+1):(2*n)]    # calculation of x[i]-y[i]
		W.max<-n*(n+1)/2        # some statistics of the null density; W.min is always 0
		W.mean<-n*(n+1)/4
		if(sum(Y==0)>0 & zero.rand){    # if zero.rand=TRUE, any zero (i.e. any x[i]-y[i]) will be set to
			cat("zeros:",sum(Y==0),"\n")   # a very small positive or negative value. Our way of randomly
			Y[which(Y==0)]<-sample(c(1e-008,-1e-008),sum(Y==0),replace=TRUE)  # assigning a sign to zeros
		}
		y.wilk<-NULL
		for(i in 1:n.genes)   # here the null density of a Wilcoxon Sign-Rank Test is calculated
			y.wilk[i]<-sum(rank(abs(Y[i,]))[Y[i,]>0])
		f.null<-dsignrank(0:W.max,n)
	}
	if(sum(y.wilk!=round(y.wilk))>0){   # are there any ties?
		cat("tied Wilcoxon scores:", sum(y.wilk!=round(y.wilk)),"\n","\n")
		if(!ties.rand){  # a conservative way is used to get rid of non-integer test statistics
			y.wilk[which(y.wilk!=round(y.wilk) & y.wilk>W.mean)]<-floor(y.wilk[which(y.wilk!=round(y.wilk) & y.wilk>W.mean)])
			y.wilk[which(y.wilk!=round(y.wilk) & y.wilk<W.mean)]<-ceiling(y.wilk[which(y.wilk!=round(y.wilk) & y.wilk<W.mean)])
		}
		if(ties.rand){ # non-integer teststatistics are randomly assigned to either the next lower or upper
			y.rand<-sample(c(-0.5,0.5),length(which(y.wilk!=round(y.wilk))),replace=TRUE)            # integer
			y.wilk[which(y.wilk!=round(y.wilk))]<-y.wilk[which(y.wilk!=round(y.wilk))]+y.rand
	}}
	W<-as.numeric(names(table(y.wilk)))

	if(length(W)!=length(f.null)){
		if(!paired)
			f.null<-dwilcox(W-W.min,n.x,n.y)
		if(paired)
			f.null<-f.null[W+1]
	}
	count<-as.numeric(table(y.wilk))   # number of observations for each possible test statistic value
	if(!ties.rand & W.mean!=round(W.mean)){  # in the conservative way of tie-breaking it could happen that
		m<-which(W==W.mean)                     # there are still ties - if W.mean is a non-integer
		tie.count<-count[m]/2                  # to get rid of these non-integer values, they will be equally
		count[m-1]<-count[m-1]+tie.count        # distributed to either the next lower or upper integer
		count[m+1]<-count[m+1]+tie.count
		W<-W[-m]
		count<-count[-m]
	}
	offset.value<-if(use.offset) log(f.null) else rep(0,length(count))
	glm.out<-glm(count~ns(W,ns.df)+offset(offset.value),family=poisson)  # a poisson regression with natural splines
																	             # and offset is used
	f.x<-glm.out$fitted/n.genes   # to estimate the density of the observed Ws

	plot(W,n.genes*f.null,type="l",xlab="Wilcoxon score",ylab="number of genes",main="Null and mixture density")
	lines(W,n.genes*glm.out$fitted,lty=4)      # a plot of the null and the mixture density
	points(W,count)       # with the observed tie-breaked counts for each possible value of W
	legend(min(W)-1,n.genes*max(f.null,f.x),bty="n",c("null","mixture","observed"),cex=.8,lty=c(1,2,0),pch=c(-1,-1,1))

	if(is.na(p0)){
		vec.lambda<-NULL
		vec.p0 <- NULL
		spline.out<-NULL
		if(stable.p0){  # p0 is estimated in a more stable way using our interpretation of the suggestion
			  # of Remark F in Efron et al.(2001)

			for(i in 1:(floor(length(f.null-1)/2)+1)){
				vec.p0[i]<-sum(f.x[i:(length(f.null)-i+1)])/sum(f.null[i:(length(f.null)-i+1)])
				vec.lambda[i]<-1-sum(f.null[i:(length(f.null)-i+1)])
			}
			weights<-if(!use.weights) rep(1,length(vec.lambda)) else 1-vec.lambda
			spline.out<-smooth.spline(vec.lambda,vec.p0,w=weights,df=3)
			p0<-min(predict(spline.out,1)$y,1)
		}
		if(!stable.p0)     # the simple p0-estimation
			p0<-min(f.x/f.null)
	}

	cat("p0:",round(p0,4),"\n")
	p1<-1-p0*f.null/f.x      # calculation of the posterior
	x11()   # plot of the posterior
	plot(W,p1,type="b",xlab="Wilcoxon score",ylab="Pr(gene significant|W score)",main="Posterior probabilities")
	abline(h=delta,lty=4)
	points(W[which(p1>=delta)],p1[which(p1>=delta)],col=3)  # mark the "significant" W scores
	tab.sig<-table(y.wilk)[which(p1>=delta)]  # table of the "significant" W scores
	nsig<-sum(tab.sig)      # number of the significant genes
	false<-sum(f.null[which(p1>=delta)])*length(y.wilk)  # number of falsely called genes
	fdr<-p0*false/max(nsig,1)  # FDR
	cat("Number of significant genes:",nsig,"\n")  # some output
	ebam.out<-NULL

	if(nsig>0){
		cat("falsely called genes:",round(false,2),"\n","FDR:",round(fdr,4),"\n")
		p1.sig<-p1[which(p1>=delta)]
		local<-1-p1.sig      # local FDR for the significant genes
		q.value<-q.value.wilc(y.wilk,p0,n.x,n.y,paired=paired)$q.value
		row.sig.genes<-NULL
		for(i in 1:length(tab.sig))
			row.sig.genes<-c(row.sig.genes,which(y.wilk==names(tab.sig)[i])) # make some output
		ebam.output<-cbind("W"=y.wilk[row.sig.genes],"q-value"=round(q.value[row.sig.genes],4))
		if(R.fold){
			fold.change<-R.fold.cal(R.dataset[row.sig.genes,],x,y,na.rm=na.rm)
			ebam.output<-cbind(ebam.output,"R-fold"=round(fold.change[,3],4))
		}
		if(!is.na(col.gene.name))
			ebam.output<-cbind(ebam.output,"gene"=substring(data[row.sig.genes,col.gene.name],1,50))
		if(!is.na(col.accession))
			ebam.output<-cbind("access"=data[row.sig.genes,col.accession],ebam.output)
		ebam.output<-cbind("ID"=row.sig.genes,ebam.output)
		ebam.out<-as.data.frame(rbind(tab.sig,round(p1.sig,4),round(local,4)))
		ebam.out<-cbind(c("genes","p1","local"),ebam.out)
		names(ebam.out)<-c(" ",names(tab.sig))

		cat("\n","Wilcoxon test statistics of significant genes:","\n")
		print(ebam.out)
		if(!is.na(file.out)){  # store some output in a file
			cat("Results of the empirical Bayes Analysis of Microarrays using Wilcoxon Rank Statistics","\n","\n",
			"\n","Significance criterion: p1 >=",delta,"\n","\n","p0:",round(p0,4),"\n","significant genes:",nsig,
			"\n","falsely called genes:",round(false,4),"\n","FDR:",round(fdr,4),"\n","\n","Wilcoxon Rank Sums of significant genes:","\n","\n",file=file.out)
			write.table(t(dimnames(ebam.out)[[2]]),file=file.out,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
			write.table(ebam.out,file=file.out,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
			cat("\n","\n","\n","Genes called significant:","\n","\n",file=file.out,append=TRUE)
			write.table(t(dimnames(ebam.output)[[2]]),file=file.out,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
			write.table(ebam.output,file=file.out,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
			cat("\n","\n","Output is stored in",file.out,"\n")
		}
	}
	mat.out<-cbind(W,count,f.null,f.x,p1)
	structure(list(nsig=nsig,false=false,fdr=fdr,ebam.out=ebam.out,mat.out=mat.out,p0=p0,
		glm.out=glm.out,f.x=f.x,f.null=f.null,vec.p0=vec.p0,vec.lambda=vec.lambda,
		y.wilk=y.wilk,spline.out=spline.out,ebam.output=ebam.output,row.sig.genes=row.sig.genes))
}
