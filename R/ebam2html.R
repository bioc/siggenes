ebam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,findA0=NULL,
		varName=NULL,entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=FALSE,ug=FALSE,
		fullname=FALSE,chipname="",cdfname=NULL,which.refseq="NM",refsnp=NULL,
		max.associated=2,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),plotFindArgs=plotFindArguments(),
		bg.plot.adjust=FALSE,plotname=NULL,plotborder=0,tableborder=1,
		new.window=TRUE,load=TRUE,...){
	if(!is(object,"EBAM"))
		stop("object must be an object of class EBAM.")
	if(!is.null(findA0) && !is(findA0,"FindA0"))
		stop("findA0 must be an object of class FindA0.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		findA0=findA0,varName=varName,entrez=entrez,refseq=refseq,symbol=symbol,
		omim=omim,ug=ug,fullname=fullname,chipname=chipname,cdfname=cdfname,
		refsnp=refsnp,max.associated=max.associated,n.digits=n.digits,bg.col=bg.col,
		text.col=text.col,link.col=link.col,plotArgs=plotArgs,plotFindArgs=plotFindArgs,
		which.refseq=which.refseq,bg.plot.adjust=bg.plot.adjust,plotname=plotname,
		plotborder=plotborder,tableborder=tableborder,new.window=new.window,load=load,...)
}

