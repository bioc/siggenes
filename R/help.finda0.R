`help.finda0` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","plot"))
		stop("'method' must be either print or plot.")
	fp<-find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"finda0.html",sep="."))
	browseURL(fp)
}

