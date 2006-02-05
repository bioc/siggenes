help.sam<-function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method %in% c("print","summary","plot","identify"))
		stop("'method' must be either print, summary, plot or identify.")
	fp<-.find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"sam.html",sep="."))
	browseURL(fp)
}


