d.null<-function(X,mat.samp,d,type.mt,s0,med=FALSE,n.subset=10){
	if(med)
		n.subset<-1
	n.samp<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,n.samp,n.subset),n.samp+1))
	n.int<-length(seq.samp)-1
	n.row<-nrow(X)
	d.mat<-mat.neg<-mat.pos<-matrix(0,n.row,n.int)
	le.cl<-ncol(mat.samp)/ifelse(type.mt=="pairt",2,1)
	if(type.mt=="pairt"){
		X<-X[,2*(1:le.cl)]-X[,2*(1:le.cl)-1]
		mat.samp<-mat.samp[,2*(1:le.cl)]
	}
	d.rank<-rank(-abs(d),ties="first")
	for(j in 1:n.int){
		tmp<-mat.samp[seq.samp[j]:(seq.samp[j+1]-1),,drop=FALSE]
		#if(!is.matrix(tmp))
		#	tmp<-matrix(tmp,1)
		dperm.out <- build.dperm(X,tmp,type.mt,s0,n.row,le.cl)	
		d.mat[,j] <- rowMeans(as.matrix(dperm.out), na.rm=TRUE) * nrow(tmp)
		mat.pos[,j] <- rank(-c(dperm.out[dperm.out>=0],abs(d)), na.last=NA,
			ties="first")[sum(dperm.out>=0, na.rm=TRUE) + (1:n.row)] - d.rank
		mat.neg[,j] <- rank(c(dperm.out[dperm.out<0],-abs(d)), na.last=NA,
			ties="first")[sum(dperm.out<0, na.rm=TRUE) + (1:n.row)] - d.rank
	}
	B<-nrow(mat.samp)
	d.bar<-rowSums(as.matrix(d.mat))/B
	p.value<-(rowSums(as.matrix(mat.pos))+rowSums(as.matrix(mat.neg)))/(n.row*B)
	vec.false<-numeric(n.row)
	if(med){
		vec.false[d>=0]<-apply(as.matrix(mat.pos[d>=0,]),1,median)
		if(type.mt!="f")
			vec.false[d<0]<-apply(as.matrix(mat.neg[d<0,]),1,median)
	}
	else{
		vec.false[d>=0]<-rowSums(as.matrix(mat.pos[d>=0,]))/B
		if(type.mt!="f")
			vec.false[d<0]<-rowSums(as.matrix(mat.neg[d<0,]))/B
	}
	invisible(return(list(d.bar=d.bar,p.value=p.value,vec.false=vec.false,d.mat=d.mat,
		mat.pos=mat.pos,mat.neg=mat.neg)))
} 
 

