find.a0<-function(data,cl,B=100,balanced=FALSE,mat.samp=NULL,delta=0.9,alpha=(0:9)/10,
		include.0=TRUE,p0=NA,plot.legend=TRUE,na.rm=FALSE,rand=TRUE){
	X.name<-match.call()$data
	if(is(data,"exprSet")){
		require(affy)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
	}
	if(is.data.frame(cl) || is.matrix(cl)) 
        	cl<-pairt.cl.transform(cl,ncol(data))
	xy.out<-xy.cal(cl,FALSE,TRUE)
	x<-xy.out$x
	y<-xy.out$y
	paired<-xy.out$paired
	tmp<-find.a0.old(data,x,y,paired=paired,mat.samp=mat.samp,B=B,balanced=balanced,na.rm=na.rm,
		delta=delta,alpha=alpha,include.0=include.0,p0=p0,rand=rand,plot.legend=plot.legend)
	invisible(structure(list(X.name=X.name,r=tmp$r,s=tmp$s,r.perm=tmp$r.perm,s.perm=tmp$s.perm,
		mat.samp=tmp$mat.samp,sig.a0=tmp$sig.a0,a0=tmp$a0,delta=tmp$delta,vec.a0=tmp$vec.a0,
		x=x,y=y)))
}
 
 

