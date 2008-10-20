sam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,varName=NULL,
		entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=FALSE,ug=FALSE,fullname=FALSE,bonf=FALSE,
		chipname="",cdfname=NULL,which.refseq="NM",refsnp=NULL,max.associated=2,
		n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),bg.plot.adjust=FALSE,plotname=NULL,
		plotborder=0,tableborder=1,new.window=TRUE,load=TRUE,...){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	if(length(delta)!=1)
		stop("delta must be a numeric value.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		varName=varName,entrez=entrez,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		fullname=fullname,chipname=chipname,cdfname=cdfname,refsnp=refsnp,bonf=bonf,
		max.associated=max.associated,n.digits=n.digits,bg.col=bg.col,text.col=text.col,
		link.col=link.col,plotArgs=plotArgs,which.refseq=which.refseq,
		bg.plot.adjust=bg.plot.adjust,plotname=plotname,plotborder=plotborder,
		tableborder=tableborder,new.window=new.window,load=load,...)
}




	