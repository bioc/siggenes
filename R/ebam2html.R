ebam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,findA0=NULL,
		varName=NULL,ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,refsnp=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),plotFindArgs=plotFindArguments(),
		bg.plot.adjust=FALSE,plotname=NULL,plotborder=0,tableborder=1,
		new.window=TRUE,...){
	if(!is(object,"EBAM"))
		stop("object must be an object of class EBAM.")
	if(!is.null(findA0) && !is(findA0,"FindA0"))
		stop("findA0 must be an object of class FindA0.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		findA0=findA0,varName=varName,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,refsnp=refsnp,n.digits=n.digits,bg.col=bg.col,
		text.col=text.col,link.col=link.col,plotArgs=plotArgs,plotFindArgs=plotFindArgs,
		bg.plot.adjust=bg.plot.adjust,plotname=plotname,plotborder=plotborder,
		tableborder=tableborder,new.window=new.window,...)
}

