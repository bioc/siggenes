`help.ebam` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot"))
		stop("'method' must be either print, summary, or plot.")
	fp<-find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"ebam.html",sep="."))
	browseURL(fp)
}

