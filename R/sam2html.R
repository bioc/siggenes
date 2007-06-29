sam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,varName=NULL,
		ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,refsnp=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),bg.plot.adjust=FALSE,plotname=NULL,
		plotborder=0,tableborder=1,new.window=TRUE,...){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	if(length(delta)!=1)
		stop("delta must be a numeric value.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		varName=varName,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,refsnp=refsnp,n.digits=n.digits,bg.col=bg.col,
		text.col=text.col,link.col=link.col,plotArgs=plotArgs,
		bg.plot.adjust=bg.plot.adjust,plotname=plotname,plotborder=plotborder,
		tableborder=tableborder,new.window=new.window,...)
}




	