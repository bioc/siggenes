`compFailureSubset` <-
function(data,mat.samp,z.mat,ints,z.norm,ints.norm,FUN,vec.a0,
		z.range,n.interval=139){
	B<-nrow(mat.samp)
	mat.r<-mat.s<-matrix(0,nrow(data),B)
	n.a0<-length(vec.a0)
	tmp.fail<-tmp.fail.norm<-matrix(0,n.interval,n.a0)
	for(i in 1:B){
		fun.out<-FUN(data,mat.samp[i,])
		mat.r[,i]<-fun.out$r
		mat.s[,i]<-fun.out$s
	}
	for(i in 1:n.a0){
		tmp<-as.vector(mat.r/(mat.s+vec.a0[i]))
		tmp.fail[,i]<-getFailure(tmp,z.mat[,i],ints[[i]],z.range=z.range[,i])
		tmp.fail.norm[,i]<-getFailure(tmp,z.mat[,i],ints.norm,z.norm=z.norm,
			n.interval=n.interval)
	}
	return(list(tmp.fail=tmp.fail,tmp.fail.norm=tmp.fail.norm))
}

