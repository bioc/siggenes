`args.ebam` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot"))
		stop("'method' must be either print, summary, or plot.")
	if(method=="print")
		cat("function (x, delta = NULL, n.digits = 4)\n")
	if(method=="summary")
		cat("function(object, delta = NULL, n.digits = 4, what = \"both\",",
			"ll = FALSE,\n    chip = \"\", file = \"\", sep = \"\\t\",",
			"quote = FALSE, dec = \".\")\n")
	if(method=="plot")
		cat("function (x, y, pos.stats = NULL, sig.col = 3, sig.cex = 1,",
			"pch = NULL,\n    cexStats = 0.8, main = NULL, xlab = NULL,",
			"ylab = NULL, y.intersp = 1.3, ...)\n")
}

