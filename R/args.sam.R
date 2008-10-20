args.sam<-function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot","identify"))
		stop("'method' must be either print, summary, plot or identify.")
	if(method=="print")
		cat("function (x, delta = NULL, n.digits = 3)\n")
	if(method=="summary")
		cat("function(object, delta = NULL, n.digits = 5, what = \"both\",",
			"entrez = FALSE,\n    bonf = FALSE, chip = \"\", file = \"\", sep = \"\\t\",",
			"quote = FALSE, dec = \".\")\n")
	if(method=="plot")
		cat("function (x, y, pos.stats = NULL, sig.col = 3, xlim = NULL,",
			"ylim = NULL,\n    main = NULL, xlab = NULL, ylab = NULL,",
			"pty = \"s\", lab = c(10, 10, 7),\n    pch = NULL, sig.cex = 1,",
			"helplines = TRUE,","...)\n")
	if(method=="identify")
		cat("function (x, showText = TRUE, getInfo = TRUE, pos = 4, cex = 0.8,\n",
			"   add.xy = numeric(2), n.digits = 4, ask = FALSE,",
			"entrez = FALSE,\n    browse = FALSE, chip= \"\", ...)\n")
}


