.mt.BLIM<-2^30
.mt.naNUM<- -93074815


d.stat<-function(data,cl,var.equal=FALSE,B=100,med=FALSE,s0=NA,s.alpha=seq(0,1,0.05),
		include.zero=TRUE,n.subset=10,mat.samp=NULL,B.more=0.1,B.max=30000,
		gene.names=NULL,R.fold=1,use.dm=FALSE,R.unlog=TRUE,na.replace=TRUE,
		na.method="mean",rand=NA){
	data<-as.matrix(data)
	mode(data)<-"numeric"
	n.genes<-nrow(data)
	excluded.genes<-logical(n.genes)
	if(any(rowSums(is.na(data))>0)){
		na.out<-na.handling(data,na.replace=na.replace,na.method=na.method)
		data<-na.out$X
		excluded.genes[na.out$NA.genes]<-TRUE
		rm(na.out)
	}
	adjust.out<-adjust.for.mt(data,cl,var.equal=var.equal)
	X<-adjust.out$X
	cl.mt<-adjust.out$cl.mt
	type.mt<-adjust.out$type.mt
	msg<-adjust.out$msg
	if(!type.mt%in%c("t","t.equalvar"))
		R.fold<-0
	if(R.fold>0){
		mat.fold<-Rfold.cal(data,cl.mt,unlog=R.unlog,R.fold=R.fold,use.dm=use.dm)
		n.fulfill<-sum(mat.fold[,2])
		if(n.fulfill==0)
			stop("None of the variables has a fold change larger than ",R.fold,".")
		if(n.fulfill<20)
			stop("Less than 20 variables have a fold change larger than ",R.fold,".")
		if(R.fold!=1)
			msg<-c(msg,paste("Number of variables having a fold change >=",R.fold,"or <=",
				round(1/R.fold,4),":",n.fulfill,"\n\n"))
		fold.out<-mat.fold[,1]
		if(n.fulfill<nrow(mat.fold)){
			fc.genes<-which(mat.fold[,2]==0)
			fold.out<-fold.out[-fc.genes]
			X<-X[-fc.genes,]
			excluded.genes[!excluded.genes][fc.genes]<-TRUE
		}
	}			
	else
		fold.out<-numeric(0)
	rm(data,adjust.out)
	mt1.out<-mt.teststat.num.denum(X,cl.mt,test=type.mt)
	r<-mt1.out$teststat.num
	s<-mt1.out$teststat.denum
	if(any(round(s,10)==0)){
		var0.genes<-which(round(s,10)==0)
		r<-r[-var0.genes]
		s<-s[-var0.genes]
		X<-X[-var0.genes,]
		if(R.fold>0)
			fold.out<-fold.out[-var0.genes]
		excluded.genes[!excluded.genes][var0.genes]<-TRUE
		warning("There are ",length(var0.genes)," variables with zero variance. These variables are removed,",
			"\n","and their d-values are set to NA.",call.=FALSE)
	}
	if(is.na(s0)){
		s0.out<-fudge2(r,s,alpha=s.alpha,include.zero=include.zero)
		s0<-s0.out$s.zero
		msg<-c(msg,s0.out$msg)
	}
	else
		msg<-c(msg,paste("s0 =",round(s0,4),"\n\n"))
	d<-r/(s+s0)
	#d.sort<-sort(d)
	mat.samp<-setup.mat.samp(cl.mt,type.mt,B=B,mat.samp=mat.samp,B.more=B.more,B.max=B.max,rand=rand)
	B.full<-round(ifelse(type.mt=="pairt",2^(length(cl.mt)/2),
		exp(lgamma(length(cl.mt)+1)-sum(lgamma(table(cl.mt)+1)))))
	msg<-c(msg,paste("Number of permutations:",nrow(mat.samp),
		ifelse(nrow(mat.samp)==B.full,"(complete permutation)",""),"\n\n"),
		paste(ifelse(med,"MEDIAN","MEAN"),"number of falsely called variables is computed.\n\n"))
	dnull.out<-d.null(X,mat.samp,d,type.mt,s0,med=med,n.subset=n.subset)
	d.bar<-dnull.out$d.bar
	if(nrow(X)==n.genes){
		p.value<-dnull.out$p.value
		vec.false<-dnull.out$vec.false
	}
	else{
		p.value<-vec.false<-d.new<-s.new<-rep(NA,n.genes)
		p.value[!excluded.genes]<-dnull.out$p.value
		vec.false[!excluded.genes]<-dnull.out$vec.false
		d.new[!excluded.genes]<-d
		s.new[!excluded.genes]<-s
		d<-d.new
		s<-s.new
		if(R.fold>0){
			f.new<-rep(NA,n.genes)
			f.new[!excluded.genes]<-fold.out
			fold.out<-f.new
		}
			
	}
	rm(dnull.out)
	if(!is.null(gene.names))
		names(d)<-names(p.value)<-names(vec.false)<-names(s)<-substring(gene.names,1,50)
	invisible(list(d=d,d.bar=d.bar,p.value=p.value,vec.false=vec.false,discrete=FALSE,s=s,s0=s0,
		mat.samp=mat.samp,msg=msg,fold=fold.out))	
}
 
 

