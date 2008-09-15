link.siggenes<-function(object,delta,filename,gene.names=NULL,addDataFrame=TRUE,
		entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=FALSE,ug=FALSE,fullname=FALSE,
		which.refseq="NM",chipname="",cdfname=NULL,n.digits=3,title=NULL,bg.col="white",
		text.col="black",link.col="blue",tableborder=1,new.window=TRUE,load=TRUE){
	if(any(c(entrez,refseq,symbol,omim,ug,fullname)))
		chipname<-check.chipname(chipname,object@chip,cdfname)
	genenames<-list.siggenes(object,delta,gene.names=gene.names)
	dataframe<-if(addDataFrame) 
			pretty.mat.sig(summary(object,delta,what="genes")@mat.sig[,-1],digits=n.digits) 
		else NULL
	link.genes(genenames,filename,entrez=entrez,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		fullname=fullname,which.refseq=which.refseq,chipname=chipname,cdfname=cdfname,
		dataframe=dataframe,title=title,bg.col=bg.col,text.col=text.col,
		link.col=link.col,tableborder=tableborder,new.window=new.window,load=load)
}
	


		