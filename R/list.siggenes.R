`list.siggenes` <-
function(object,delta,file="",gene.names=NULL,order=TRUE,text=NULL,append=FALSE){
	if(!class(object)%in%c("SAM","EBAM"))
		stop("object must be either a SAM or an EBAM object.")
	if(length(delta)!=1)
		stop("delta must have length 1.")
	if(delta<=0)
		stop("delta must be larger than 0.")
	if(is(object,"EBAM") & delta>1)
		stop("delta must be smaller than 1.")
	d<-if(is(object,"SAM")) object@d else object@z
	if(is.null(gene.names))
		gene.names<-names(d)
	if(is.null(gene.names))
		stop("gene.names must be specified.")
	if(length(gene.names)!=length(d))
		stop("The length of gene.names differs from the number of genes.")
	if(!is.null(names(d)) & any(gene.names!=names(d)))
		stop("Some of the gene.names differ from the gene names of the ",
			ifelse(is(object,"SAM"),"SAM","EBAM")," object.")
	rsg<-summary(object,delta=delta,what="stats")@row.sig.genes
	if(order)
		rsg<-rsg[rev(order(abs(d[rsg])))]
	sig.genes<-gene.names[rsg]
	if(file!=""){
		if(is.null(text))
			text<-paste(object@msg[1],"Delta = ",delta,"\n\n",sep="")
		if(text=="")
			text<-NULL
		cat(text,if(!is.null(text)) "\n",file=file,append=append,sep="")
		write(sig.genes,file=file,append=TRUE)
	}
	else
		sig.genes
}

