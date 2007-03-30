`compFailureMat2` <-
function(data,mat.samp,z.mat,ints,z.norm,ints.norm,FUN,vec.a0,n.chunk=1,
		n.interval=139){
	B<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,B,ceiling(B/n.chunk)),B+1))
	n.seq<-length(seq.samp)-1
	mat.failure<-mat.fail.norm<-matrix(0,n.interval,length(vec.a0))
	z.range<-apply(z.mat,2,range)
	for(j in 1:n.seq){
		tmp.samp<-mat.samp[seq.samp[j]:(seq.samp[j+1]-1),]
		tmp<-compFailureSubset(data,tmp.samp,z.mat,ints,z.norm,ints.norm,FUN,
			vec.a0,z.range,n.interval=n.interval)
		mat.failure<-mat.failure+tmp$tmp.fail
		mat.fail.norm<-mat.fail.norm+tmp$tmp.fail.norm
	}
	colnames(mat.failure)<-colnames(mat.fail.norm)<-round(vec.a0,4)
	return(list(mat.failure=mat.failure,mat.fail.norm=mat.fail.norm))
}

