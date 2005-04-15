link.siggenes<-function(object,delta,filename,gene.names=NULL,addDataFrame=TRUE,
		ll=FALSE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,title=NULL,bg.col="white",text.col="black",
		link.col="blue",tableborder=1,new.window=TRUE){
	if(any(c(ll,refseq,symbol,omim,ug)))
		chipname<-check.chipname(chipname,object@chip,cdfname)
	genenames<-list.siggenes(object,delta,gene.names=gene.names)
	dataframe<-if(addDataFrame) round(summary(object,delta,what="genes")@mat.sig[,-1],
		n.digits) else NULL
	link.genes(genenames,filename,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,dataframe=dataframe,
		title=title,bg.col=bg.col,text.col=text.col,link.col=link.col,
		tableborder=tableborder,new.window=new.window)
}
	


		