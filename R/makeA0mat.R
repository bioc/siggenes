`makeA0mat` <-
function(z.norm,mat.post,vec.p0,vec.a0,B,delta=0.9){
	n.a0<-length(vec.a0)
	mat<-matrix(0,n.a0,2)
	colnames(mat)<-c("Number","FDR")
	for(i in 1:n.a0)
		mat[i,]<-compNumber(z.norm,mat.post[,i],vec.p0[i],B,
			delta=delta)[,2:3]
	out<-data.frame(a0=round(vec.a0,4),Quantile=names(vec.a0),round(mat,4))
	rownames(out)<-1:nrow(out)
	ids<-which.max(mat[,1])
	list(tab=out,suggest=vec.a0[ids])
}

