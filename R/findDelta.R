`findDelta` <-
function(object,fdr=NULL,genes=NULL,prec=6,initial=NULL,verbose=FALSE){
	if(!class(object)[1] %in% c("SAM","EBAM"))
		stop("object must be an object of class SAM or EBAM.")
	isSAM<-is(object,"SAM")[1]
	if(is.null(fdr) & is.null(genes))
		stop("Either fdr or genes must be specified.")
	if(!is.null(fdr) & !is.null(genes))
		stop("Only one of fdr and genes can be specified.")
	delta<-checkInitialDelta(object,initial=initial,isSAM=isSAM,prec=prec)
	if(!is.null(fdr)){
		if(fdr<=0 | fdr>=1)
			stop("fdr must be larger than 0 and smaller than 1.")
		out<-findNumber(object,fdr,delta=delta,isSAM=isSAM,prec=prec,verbose=verbose)
	}
	else{
		if(genes!=round(genes))
			stop("genes must be an integer.")
		tmp<-ifelse(isSAM,"d","z")
		n.genes<-length(slot(object,tmp))
		if(genes<=0 | genes>n.genes)
			stop("genes must be an integer between 0 and ",n.genes,".")
		out<-findFDR(object,genes,delta=delta,isSAM=isSAM,prec=prec,verbose=verbose)
	}
	invisible(out)
}

