# Copyright (C) 2003 Holger Schwender

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


sam.wilc<-function(data,cl,delta=1:max(abs(W.diff)),na.rm=FALSE,zero.rand=TRUE,
		rand=NA,graphic=TRUE,thres=round(quantile(2:max(abs(W.diff)),(0:3)/3)),
		use.numbers=TRUE){
	X.name<-match.call()$data
	xy.out<-xy.cal(cl,TRUE)
	x<-xy.out$x
	y<-xy.out$y
	paired<-xy.out$paired
	wilc.out<-wilc.cal(data,x,y,paired=paired,zero.rand=zero.rand,rand=rand,na.rm=na.rm)
   	W<-wilc.out$W        # observed unsorted W-values
    	W.sort<-sort(W)      # observed sorted W-values
    	W.exp<-wilc.out$W.exp   # expected (sorted) W-values
    	var.0.genes<-wilc.out$var.0.genes  # index of genes with variance Zero
    	n.genes<-sum(!is.na(W))    # number of genes with non missing W-value
    	W.exp.number<-n.genes*wilc.out$f.null    # vector of expected numbers of the W-values under the null
    	W.exp.value<-as.numeric(names(W.exp.number))    # vector of the possible W-values
    	n.exp<-length(W.exp.number)  
    	vec.lambda<-NULL
    	vec.p0<-NULL
    	for(i in 1:floor(n.exp/2)){   # estimation of p0
        	vec.p0[i]<-sum(W>=W.exp.value[i] & W<=W.exp.value[n.exp-i+1],na.rm=TRUE)/sum(W.exp.number[i:(n.exp-i+1)])
        	vec.lambda[i]<-1-sum(W.exp.number[i:(n.exp-i+1)])/n.genes
    	}
    	library(modreg)
    	spline.out<-smooth.spline(vec.lambda,vec.p0,w=1-vec.lambda,df=3)
    	p0<-min(predict(spline.out,1)$y,1)    
    	W.diff<-W.sort-W.exp   # computation of W.diff  (we are looking for |W.diff|>=delta)
    	if(any(delta<=0) || any(delta!=round(delta)))
		stop("Delta must be a positive integer.")
	table.W.exp<-cbind(W.exp.value,W.exp.number)
    	sam.fdr.out<-sam.fdr(W.sort,W.exp,table.W.exp,p0,delta=delta,wilc=TRUE,graphic=FALSE)
    	FDR<-sam.fdr.out$mat.fdr
    	print(sam.fdr.out$tab.fdr)
    
    	table.count<-table(W.exp,W.sort)  # for the SAM Plot the numbers of observations which correspond to the
    	W.exp.value<-as.numeric(dimnames(table.count)[[1]])  # possible points is computed
    	W.value<-as.numeric(dimnames(table.count)[[2]])
    	n.exp<-length(W.exp.value)
    	mat.count<-matrix(c(rep(W.exp.value,length(W.value)),rep(W.value,each=n.exp), 
        	as.vector(table.count)),ncol = 3)
    	mat.count<-mat.count[-which(mat.count[,3]==0),]
    	n.sig<-NULL    
    	# SAM Plot
    	if(graphic){
        	sam.plotter(mat.count[,2],mat.count[,1],thres,main="SAM Plot using Wilcoxon Rank Statistics",
            		color=2:(length(thres)+1),make.legend=TRUE,wilc=TRUE,use.numbers=use.numbers,count=mat.count[,3])
        	X11() 
        	# delta vs. FDR and delta vs. #significant genes
        	roller.coaster(FDR) # Delta vs. FDR and Delta vs. #significant genes
	}
    	invisible(structure(list(X.name=X.name,W=W,W.exp=W.exp,W.sort=W.sort,W.exp.number=W.exp.number,
		p0=p0,spline.out=spline.out,FDR=FDR,mat.count=mat.count,
		x=x,y=y,paired=paired,use.numbers=use.numbers,rand=rand)))
}