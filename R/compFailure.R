`compFailure` <-
function(data,mat.samp,z,interval,a0=0,type.mt=NULL,n.subset=5,fast=FALSE,n.cat=NULL){
	B<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,B,ceiling(B/n.subset)),B+1))
	n.int<-length(interval)-1
	n.seq<-length(seq.samp)-1
	vec.fail<-numeric(n.int)
	n.row<-nrow(data)
	le.cl<-ncol(mat.samp)
	z.range<-range(z)
	if(!fast){
		vec.pos<-vec.neg<-numeric(n.row)
		z.rank<-rank(-abs(z),ties="first")
	}
	else
		vec.pos<-vec.neg<-NULL
	for(i in 1:n.seq){
		tmp.samp<-mat.samp[seq.samp[i]:(seq.samp[i+1]-1),]
		if(is.null(n.cat))
			z.perm<-build.dperm(data,tmp.samp,type.mt,a0,n.row,le.cl)
		else
			z.perm<-compPermStat(data,tmp.samp,n.cat)
		z.perm<-as.vector(z.perm)
		tmp<-getFailure(z.perm,z,interval,z.range=z.range,n.interval=n.int)
		vec.fail<-vec.fail+tmp
		if(!fast){
			tmp<-compFalse(z,z.perm,z.rank,n.row)
			vec.pos<-vec.pos+tmp$vec.pos
			vec.neg<-vec.neg+tmp$vec.neg
		}
	}
	return(list(vec.fail=vec.fail,vec.pos=vec.pos,vec.neg=vec.neg))
}

