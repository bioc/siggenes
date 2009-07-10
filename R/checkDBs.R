checkDBs <- function(entrez, refseq, symbol, omim, ug, fullname, chipname, load=TRUE){
	getAnnMap("", chipname, load=load)
	obj <- search()
	pkg <- obj[obj %in% paste("package:", chipname, c(".db", ""), sep="")]
	vecDBs <- c("ENTREZID", "REFSEQ", "SYMBOL", "OMIM", "UNIGENE", "GENENAME")
	vecExist <- paste(chipname, vecDBs, sep="") %in% ls(pkg[1])
	pkg <- gsub("package:", "", pkg[1])
	if(entrez && !vecExist[1]){
		warning("ENTREZID does not seem to be available in ",  pkg, ".", call.=FALSE)
		entrez <- FALSE
	}
	if(refseq && !vecExist[2]){
		warning("REFSEQ does not seem to be available in ",  pkg, ".", call.=FALSE)
		refseq <- FALSE
	}
	if(symbol && !vecExist[3]){
		warning("SYMBOL does not seem to be available in ",  pkg, ".", call.=FALSE)
		symbol <- FALSE
	}
	if(omim && !vecExist[4]){
		warning("OMIM does not seem to be available in ",  pkg, ".", call.=FALSE)
		omim <- FALSE
	}
	if(ug && !vecExist[5]){
		warning("UNIGENE does not seem to be available in ",  pkg, ".", call.=FALSE)
		ug <- FALSE
	}
	if(fullname && !vecExist[6]){
		warning("GENENAME does not seem to be available in ",  pkg, ".", call.=FALSE)
		fullname <- FALSE
	}
	return(list(entrez=entrez, refseq=refseq, symbol=symbol, omim=omim, ug=ug, 
		fullname=fullname))
}

