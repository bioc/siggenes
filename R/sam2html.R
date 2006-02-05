sam2html<-function(sam.out,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,
		ll=FALSE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),bg.plot.adjust=FALSE,plotname=NULL,
		plotborder=0,tableborder=1,new.window=TRUE,...){
	if(!is(sam.out,"SAM"))
		stop("'sam.out' must be an object of class SAM.")
	if(length(delta)!=1)
		stop("'delta' must be a numerical value.")
	if(any(c(ll,refseq,symbol,omim,ug))){
		if(is.null(names(sam.out@d))){
			ll<-refseq<-symbol<-omim<-ug<-FALSE
			warning("Since no gene names are specified by 'sam.out'",
				" 'll', 'refseq', 'symbol',\n","'omim' and 'ug' are set",
				" to FALSE.",call.=FALSE)
		}
		if(chipname=="" & sam.out@chip=="" & is.null(cdfname)){
			ll<-refseq<-symbol<-omim<-ug<-FALSE
			warning("Since the chip type has been specified neither",
				" by the SAM object\n","nor by 'chipname' or",
				" 'cdfname', 'll', 'refseq', 'symbol',\n",
				"'omim' and 'ug' are set to FALSE.",call.=FALSE)
		}
	}
	if(any(c(ll,refseq,symbol,omim,ug)))
		chipname<-check.chipname(chipname,sam.out@chip,cdfname)
	suffix<-tolower(substring(filename,nchar(filename)-4,nchar(filename)))
	if(suffix!=".html"){
		filename<-paste(filename,"html",sep=".")
		warning("Since the suffix of 'filename' is not 'html' '.html' is added",
			" to 'filename'.",call.=FALSE)
	}
	if(is.null(plotname))
		plotname<-paste("samplot_for_",gsub(".html","",basename(filename)),
			".png",sep="")
	file.sep<-.Platform$file.sep
	tmp<-unlist(strsplit(filename,file.sep))
	path<-if(length(tmp)>1) paste(tmp[-length(tmp)],collapse=file.sep)
		else "."	
	bg.col<-col2hex(bg.col)
	text.col<-col2hex(text.col)
	link.col<-col2hex(link.col)
	sum.out<-summary(sam.out,delta,n.digits=n.digits,chip=chipname)
	mat.fdr<-pretty.mat.fdr(sum.out@mat.fdr,digits=n.digits)
	mat.sig<-sum.out@mat.sig
	if(!is.null(mat.sig))
		mat.sig<-pretty.mat.sig(mat.sig,digits=n.digits)
	else{
		if(addGenes){
			addGenes<-FALSE
			warning("Since there are no significant genes 'addGenes' is set",
				" to FALSE.",call.=FALSE)
		}
	}
	msg<-sam.out@msg
	h2<-unlist(strsplit(msg[1],"\n"))[1]
	h2<-substring(h2,13,nchar(h2))
	outfile<-file(filename,"w")
	type<-"text/css"
	cat("<html>","<head>","<title>SAM Analysis</title>","</head>",
		paste("<body bgcolor=",bg.col," text=",text.col," link=",link.col,
		">",sep=""),"<h1 align=center> SAM Analysis </h1>",
		paste("<style type=",type,">",sep=""),
		"h3{ margin-top: 5px; margin-bottom: 2px; padding-left: 10px; text-indent: 10px;
			word-spacing: 1px}", 
		"li{ margin-top: 10px; padding-left: 0px; text-indent: 5px; word-spacing: 1px}",
		"ul{ padding-left: 50px}",	
		"</style>",
		paste("<h2 align=center>",h2,"</h2>"),
		sep="\n",file=outfile)
	if(addStats){
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<h3>General Information</h3>",sep="\n",file=outfile)
		cat("<ul>",sep="\n",file=outfile)
		if(length(sam.out@s0)==1){
			s0.msg<-msg[substring(msg,1,2)=="s0"]
			if(length(s0.msg==1))
				cat("<li>Fudge Factor: ",s0.msg,"\n",sep="",file=outfile)
		}
		cat(paste("<li>Prior Probability: p0 = ",mat.fdr[,"p0"],sep=""),
			paste("<li>Statistics for Delta = ",mat.fdr[,"Delta"],":",sep=""),
			"<ul type=\"disc\">",
			paste("<li>Number of Significant Genes:",mat.fdr[,"Called"]),
			paste("<li>Number of Falsely Called Genes:",mat.fdr[,"False"]),
			paste("<li>False Discovery Rate: ",mat.fdr[,"FDR"],sep=""),
			"</ul>","</ul>",sep="\n",file=outfile)
	}
	if(addPlot){
		suf.plot<-unlist(strsplit(plotname,"\\."))
		suf.plot<-suf.plot[length(suf.plot)]
		if(!suf.plot%in%c("jpeg","png"))
			stop("'plotname' must be either a png or a jpeg file.")
		tmp.local<-length(unlist(strsplit(plotname,file.sep)))
		if(tmp.local==1)
			plotname<-paste(path,plotname,sep=file.sep)
		FUN<-match.fun(suf.plot)
		FUN(plotname)
		if(bg.plot.adjust)
			par(bg=bg.col)
		else
			par(bg="white")
		plot(sam.out,delta,pos.stats=plotArgs$pos.stats,sig.col=plotArgs$sig.col,
			xlim=plotArgs$xlim,ylim=plotArgs$ylim,main=plotArgs$main,
			xlab=plotArgs$xlab,ylab=plotArgs$ylab,pty=plotArgs$pty,
			lab=plotArgs$lab,pch=plotArgs$pch,sig.cex=plotArgs$sig.cex,...)
		dev.off()
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<div style=\"text-align: center\"><img src=",
			if(tmp.local!=1 | path!=".") "file:///",plotname," border=",
			plotborder,"></div>","\n",sep="",file=outfile)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
	}
	if(addGenes){
		has.nonames<-all(rownames(mat.sig)==as.character(1:nrow(mat.sig)))
		if(has.nonames)
			tr<-make.tablecode(as.character(mat.sig[,"Row"]),ll=FALSE,
				refseq=FALSE,symbol=FALSE,omim=FALSE,ug=FALSE,
				chipname=chipname,dataframe=mat.sig[,-1],
				new.window=new.window,tableborder=tableborder)
		else
			tr<-make.tablecode(rownames(mat.sig),ll=ll,refseq=refseq,
				symbol=symbol,omim=omim,ug=ug,chipname=chipname,
				cdfname=cdfname,dataframe=mat.sig[,-1],
				new.window=new.window,tableborder=tableborder)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<h3 align=center>Genes Called Differentially Expressed (Using Delta = ",
			mat.fdr[,"Delta"],")</h3>","\n",sep="",file=outfile)
		cat("<p><font color=",bg.col," size=1> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<style type=text/css>",
			"p{ margin-top: 2px; margin-bottom: 2px; word-spacing: 1px}",
			"</style>",tr,sep="\n",file=outfile)
		if(has.nonames)
			cat("<p><font color=",bg.col," size=1> HALLO</font></p>",
				"<p align=center><b>Annotation:</b> Since no gene names 
				have been specified the row numbers of the genes are used 
				as names.",sep="\n",file=outfile)
	}
	cat("<p><font color=",bg.col," size=5> HALLO</font></p>","\n",sep="",
			file=outfile)	
	cat("</body>","</html>",sep="\n",file=outfile)
	close(outfile)
	cat("Output is stored in ",filename,".\n",sep="")
	if(addPlot)
		cat("The SAM Plot required by the html file is stored in ",plotname,
		".\n",sep="")
} 



	