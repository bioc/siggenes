list.siggenes<-function(object,delta,file="",gene.names=NULL,order=TRUE,text=NULL,
		append=FALSE){
	if(length(delta)!=1)
		stop("delta must have length 1.")
	if(is.null(gene.names))
		gene.names<-names(object@d)
	if(is.null(gene.names))
		stop("'gene.names' must be specified.")
	if(length(gene.names)!=length(object@d))
		stop("The length of gene.names must be equal to the number of genes.")
	if(!is.null(names(object@d)) & any(gene.names!=names(object@d)))
		stop("'gene.names' differ from the genes names of the SAM object.")
	mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
		delta=delta)
	row.sig.genes<-which(object@d>=mat.fdr[,"cutup"] | object@d<=mat.fdr[,"cutlow"])
	if(order)
		row.sig.genes<-row.sig.genes[rev(order(abs(object@d[row.sig.genes])))]
	sig.genes<-gene.names[row.sig.genes]
	if(file!=""){
		cat(text,if(!is.null(text)) "\n",file=file,append=append)
		write(sig.genes,file=file,append=TRUE)
	}
	else
		sig.genes
}
	
	
	 
 

