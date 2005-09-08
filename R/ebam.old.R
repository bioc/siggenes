ebam.old<-function(a0.out,data,a0=NA,p0=NA,delta=NA,stable=TRUE,number.int=139,local.bin=.1,
		col.accession=NA,col.gene.name=NA,q.values=TRUE,R.fold=TRUE,R.dataset=data,
		na.rm=FALSE,file.out=NA){
	r<-a0.out$r
	s<-a0.out$s
	r.perm<-a0.out$r.perm
	s.perm<-a0.out$s.perm
	paired<-a0.out$paired
	if(is.na(delta))     
		delta<-a0.out$delta
	if(is.na(a0))
		a0<-a0.out$a0
	if(a0<0)     
		stop("a0 must be larger or equal to 0.")
	if(delta>=1 || delta<=0)  
		stop("delta must be between 0 and 1.")
	B<-length(r.perm)/length(r) 
	Z.unsorted<-r/(s+a0)        
	Z<-sort(Z.unsorted)
	z.unsorted<-as.vector(r.perm/(s.perm+a0))
	z<-sort(z.unsorted)
	z[z<min(Z)]<-min(Z)    
	z[z>max(Z)]<-max(Z)    
	ratio.out<-ratio.est(Z,z,p0=p0,stable=stable,number.int=number.int)  
	mat.post<-ratio.out$mat.post
	mat.repeat<-ratio.out$mat.repeat
	p0<-ratio.out$p0
	optim.out<-ratio.out$optim.out
	sig.center<-which(mat.post[,"posterior"]>=delta) 
	nsig<-sum(mat.post[sig.center,"success"])       
	false<-sum(mat.repeat[sig.center,"n"]-mat.repeat[sig.center,"success"])/B 
	fdr<-p0*false/max(nsig,1)  
	posterior<-rep(mat.post[,"posterior"],mat.post[,"success"])  
	mat.post.Z<-cbind(Z,posterior)                               
	sig.genes<-which(mat.post.Z[,"posterior"]>=delta)    
	index <- (1:length(Z))
    	row.sig.genes <- index[order(Z.unsorted)][sig.genes]
	FDR<-c(p0=p0,nsig=nsig,false=false,fdr=fdr)   
	local2<-NULL            
	for(i in 1:length(sig.genes)){
		local2[i]<-p0*sum(z>=Z[sig.genes[i]]-local.bin & z<=Z[sig.genes[i]]+local.bin)/
			(B*sum(Z>=Z[sig.genes[i]]-local.bin & Z<=Z[sig.genes[i]]+local.bin))
	}
	mat.ebam<-cbind("Z"=round(mat.post.Z[sig.genes,1],4),"p1(Z)"=round(mat.post.Z[sig.genes,2],4),
		"local1"=round(1-mat.post.Z[sig.genes,2],4),"local2"=round(local2,4))
	if(!is.na(col.accession))
		mat.ebam<-cbind("access"=data[row.sig.genes,col.accession],mat.ebam)
	if(q.values){
		mat.qvalue<-q.value.old(Z.unsorted,z.unsorted,p0)$mat.qvalue
		q.value<-mat.qvalue[order(mat.qvalue[,1]),3]
		mat.ebam<-cbind(mat.ebam,"q-value"=round(q.value[sig.genes],5))
	}
	if(R.fold){
		fold.change<-R.fold.old(R.dataset[row.sig.genes,],a0.out$x,a0.out$y,na.rm=na.rm)
		mat.ebam<-cbind(mat.ebam,"R-fold"=round(fold.change[,3],4))
	}
	if(!is.na(col.gene.name))
		mat.ebam<-cbind(mat.ebam,"gene"=substring(data[row.sig.genes,col.gene.name],1,50))
	mat.ebam<-cbind("ID"=row.sig.genes,mat.ebam)
	ebam.out<-as.data.frame(mat.ebam)    
	cat("Using a0 =",round(a0,4),"and the original Z values, there are",nsig,
		"significant genes and",false,"falsely called genes.","\n","For p0 =",round(p0,4),
		",the FDR is",round(fdr,4),".","\n","\n")
	if(!is.na(file.out)){  
		cat("Results of the Empirical Bayes Analysis of Microarray Experiments","\n","\n","\n",
			"Significance criterion: p1(Z) >=",delta,"\n","\n","a0:",round(a0,4),"\n",
			"p0:",round(p0,4),"\n","significant genes:",nsig,"\n","falsely called genes:",
			round(false,4),"\n","FDR:",round(fdr,4),"\n","\n","\n",
			"Genes called significant:","\n","\n",file=file.out)
		write.table(t(dimnames(ebam.out)[[2]]),file=file.out,sep="     ",append=TRUE,
			row.names=FALSE,col.names=FALSE,quote=FALSE)
		write.table(ebam.out,file=file.out,sep="  ",append=TRUE,row.names=FALSE,
			col.names=FALSE,quote=FALSE)
		cat("\n","\n","local1: estimation of local FDR using the logistic regression estimates",
			"\n","local2: local FDR is estimated by binning the observations into small intervals",
			file=file.out,append=TRUE)
		cat("Output is stored in",file.out,"\n")
	}
	plot(mat.post.Z[,"Z"],mat.post.Z[,"posterior"],main="Posterior probability for Z values",
		xlab="Z values",ylab="p1(Z)")    
	abline(h=delta,lty=2)
	points(mat.post.Z[which(mat.post.Z[,2]>=delta),1],mat.post.Z[which(mat.post.Z[,2]>=delta),2],
		col=3)
	mat.Z.unsorted<-cbind(Z.unsorted,mat.post.Z[rank(Z.unsorted),2])
	invisible(structure(list(mat.repeat=mat.repeat,optim.out=optim.out,mat.post.Z=mat.post.Z,
		ebam.out=ebam.out,FDR=FDR,a0=a0,mat.Z.unsorted=mat.Z.unsorted,Z.unsorted=Z.unsorted,
		row.sig.genes=row.sig.genes,p0=p0)))
} 
 

