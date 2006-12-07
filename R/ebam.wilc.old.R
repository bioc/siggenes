ebam.wilc.old<-function(data,x,y,paired=FALSE,delta=.9,p0=NA,stable.p0=TRUE,use.offset=TRUE,
		use.weights=TRUE,ties.rand=TRUE,zero.rand=TRUE,ns.df=5,col.accession=NA,
		col.gene.name=NA,R.fold=TRUE,R.dataset=data,file.out=NA,rand=NA,na.rm=FALSE){
	require(splines)  
	if(!is.na(rand))   
		set.seed(rand)
	Y<-as.matrix(data[,c(x,y)])   
	mode(Y)<-"numeric"
	n.genes<-nrow(Y)  
	if(!paired){     
		n.x<-length(x)     
		n.y<-length(y)
		W.mean<-n.x*(n.x+n.y+1)/2    
		W.min<-n.x*(n.x+1)/2
		W.max<-n.x*(2*n.y+n.x+1)/2
		Y.rank<-t(apply(Y,1,rank))         
		y.wilk<-apply(Y.rank[,1:n.x],1,sum)
		f.null<-dwilcox(0:(n.x*n.y),n.x,n.y) 
	}
	if(paired){  
		n<-length(x)
		Y<-Y[,1:n]-Y[,(n+1):(2*n)]   
		W.max<-n*(n+1)/2       
		W.mean<-n*(n+1)/4
		if(sum(Y==0)>0 & zero.rand){   
			cat("zeros:",sum(Y==0),"\n")  
			Y[which(Y==0)]<-sample(c(1e-008,-1e-008),sum(Y==0),replace=TRUE) 		
		}
		y.wilk<-NULL
		for(i in 1:n.genes)  
			y.wilk[i]<-sum(rank(abs(Y[i,]))[Y[i,]>0])
		f.null<-dsignrank(0:W.max,n)
	}
	if(sum(y.wilk!=round(y.wilk))>0){   
		cat("tied Wilcoxon scores:", sum(y.wilk!=round(y.wilk)),"\n","\n")
		if(!ties.rand){ 
			y.wilk[which(y.wilk!=round(y.wilk) & y.wilk>W.mean)]<-floor(y.wilk[which(y.wilk!=round(y.wilk) & y.wilk>W.mean)])
			y.wilk[which(y.wilk!=round(y.wilk) & y.wilk<W.mean)]<-ceiling(y.wilk[which(y.wilk!=round(y.wilk) & y.wilk<W.mean)])
		}
		if(ties.rand){
			y.rand<-sample(c(-0.5,0.5),length(which(y.wilk!=round(y.wilk))),replace=TRUE)        
			y.wilk[which(y.wilk!=round(y.wilk))]<-y.wilk[which(y.wilk!=round(y.wilk))]+y.rand
		}
	}
	W<-as.numeric(names(table(y.wilk)))
	if(length(W)!=length(f.null)){
		if(!paired)
			f.null<-dwilcox(W-W.min,n.x,n.y)
		else
			f.null<-f.null[W+1]
	}
	count<-as.numeric(table(y.wilk))   
	if(!ties.rand & W.mean!=round(W.mean)){ 
		m<-which(W==W.mean)             
		tie.count<-count[m]/2           
		count[m-1]<-count[m-1]+tie.count 
		count[m+1]<-count[m+1]+tie.count
		W<-W[-m]
		count<-count[-m]
	}
	offset.value<-if(use.offset) log(f.null) else rep(0,length(count))
	glm.out<-glm(count~ns(W,ns.df)+offset(offset.value),family=poisson) 
	f.x<-glm.out$fitted/n.genes  
	if(is.na(p0)){
		vec.lambda<-NULL
		vec.p0 <- NULL
		spline.out<-NULL
		if(stable.p0){  			
			for(i in 1:(floor(length(f.null-1)/2)+1)){
				vec.p0[i]<-sum(f.x[i:(length(f.null)-i+1)])/sum(f.null[i:(length(f.null)-i+1)])
				vec.lambda[i]<-1-sum(f.null[i:(length(f.null)-i+1)])
			}
			weights<-if(!use.weights) rep(1,length(vec.lambda)) else 1-vec.lambda
			spline.out<-smooth.spline(vec.lambda,vec.p0,w=weights,df=3)
			p0<-min(predict(spline.out,1)$y,1)
		}
		else 
			p0<-min(f.x/f.null)   
	}
	cat("p0:",round(p0,4),"\n")
	p1<-1-p0*f.null/f.x   
	plot(W,p1,type="b",xlab="Wilcoxon score",ylab="Pr(gene significant|W score)",
		main="Posterior probabilities")
	abline(h=delta,lty=4)
	points(W[which(p1>=delta)],p1[which(p1>=delta)],col=3)  
	tab.sig<-table(y.wilk)[which(p1>=delta)]  
	nsig<-sum(tab.sig)      
	false<-sum(f.null[which(p1>=delta)])*length(y.wilk)  
	fdr<-p0*false/max(nsig,1)  
	cat("Number of significant genes:",nsig,"\n") 
	ebam.out<-NULL
	if(nsig>0){
		cat("falsely called genes:",round(false,2),"\n","FDR:",round(fdr,4),"\n")
		p1.sig<-p1[which(p1>=delta)]
		local<-1-p1.sig      
		q.value<-q.value.wilc(y.wilk,p0,n.x,n.y,paired=paired)$q.value
		row.sig.genes<-NULL 
		for(i in 1:length(tab.sig))
			row.sig.genes<-c(row.sig.genes,which(y.wilk==names(tab.sig)[i])) 
		ebam.output<-cbind("W"=y.wilk[row.sig.genes],"q-value"=round(q.value[row.sig.genes],4))
		if(R.fold){
			fold.change<-R.fold.old(R.dataset[row.sig.genes,],x,y,na.rm=na.rm)
			ebam.output<-cbind(ebam.output,"R-fold"=round(fold.change[,3],4))
		}
		if(!is.na(col.gene.name))
			ebam.output<-cbind(ebam.output,"gene"=substring(data[row.sig.genes,col.gene.name],1,50))
		if(!is.na(col.accession))
			ebam.output<-cbind("access"=data[row.sig.genes,col.accession],ebam.output)
		ebam.output<-cbind("ID"=row.sig.genes,ebam.output)	
		ebam.out<-as.data.frame(cbind(names(tab.sig),tab.sig,round(p1.sig,4),round(local,4)))
		dimnames(ebam.out)<-list(1:nrow(ebam.out),c("W","genes","p1","local"))
		if(!is.na(file.out)){  # store some output in a file
			cat("Results of the empirical Bayes Analysis of Microarrays using Wilcoxon Rank Statistics",
				"\n","\n","\n","Significance criterion: p1 >=",delta,"\n","\n","p0:",
				round(p0,4),"\n","significant genes:",nsig,"\n","falsely called genes:",
				round(false,4),"\n","FDR:",round(fdr,4),"\n","\n",
				"Wilcoxon Rank Sums of significant genes:","\n","\n",file=file.out)
			write.table(t(dimnames(ebam.out)[[2]]),file=file.out,sep=" ",append=TRUE,
				col.names=FALSE,row.names=FALSE,quote=FALSE)
			write.table(ebam.out,file=file.out,sep="      ",append=TRUE,col.names=FALSE,
				row.names=FALSE,quote=FALSE)
			cat("\n","\n","\n","Genes called significant:","\n","\n",file=file.out,
				append=TRUE)
			write.table(t(dimnames(ebam.output)[[2]]),file=file.out,sep=" ",append=TRUE,
				row.names=FALSE,col.names=FALSE,quote=FALSE)
			write.table(ebam.output,file=file.out,sep="   ",append=TRUE,row.names=FALSE,
				col.names=FALSE,quote=FALSE)
			cat("\n","\n","Output is stored in",file.out,"\n")
		}
	}
	mat.out<-cbind(W,count,f.null,f.x,p1)
	structure(list(nsig=nsig,false=false,fdr=fdr,ebam.out=ebam.out,mat.out=mat.out,p0=p0,
		glm.out=glm.out,f.x=f.x,f.null=f.null,vec.p0=vec.p0,vec.lambda=vec.lambda,
		y.wilk=y.wilk,spline.out=spline.out,ebam.output=ebam.output,row.sig.genes=row.sig.genes))
}
 
 

