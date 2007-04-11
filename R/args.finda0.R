`args.finda0` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","plot"))
		stop("'method' must be either print or plot.")
	if(method=="print")
		cat("function (x, delta = NULL)\n")
	if(method=="plot")
		cat("function (x, y, logit = TRUE, pos.legend = NULL, cexLegend = 0.8,",
			"col = NULL,\n    main = NULL, xlab = NULL, ylab = NULL,",
			"only.a0 = FALSE, lty = 1, lwd = 1,\n    y.intersp = 1.1, ...)\n")
}

