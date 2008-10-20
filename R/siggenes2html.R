siggenes2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,
		findA0=NULL,varName=NULL,entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=FALSE,ug=FALSE,
		fullname=FALSE,chipname="",cdfname=NULL,refsnp=NULL,max.associated=2,bonf=FALSE,
		n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),plotFindArgs=plotFindArguments(),
		bg.plot.adjust=FALSE,plotname=NULL,plotborder=0,tableborder=1,
		new.window=TRUE,which.refseq="NM",load=TRUE,...){
	type<-class(object)
	isSAM<-type=="SAM"
	if(length(delta)!=1)
		stop("delta must be a numeric value.")
	if(delta<=0)
		stop("delta must be larger than 0.")
	if(!isSAM && delta>=1)
		stop("delta must be smaller than 1.")
	if(any(c(entrez,refseq,symbol,omim,ug,fullname))){
		tmp<-ifelse(!isSAM,"z","d")
		if(is.null(names(slot(object,tmp)))){
			entrez<-refseq<-symbol<-omim<-ug<-fullname<-FALSE
			warning("Since no gene names are specified by 'object'",
				" 'entrez', 'refseq', 'symbol',\n","'omim', 'ug' and 'fullname' ",
				"are set to FALSE.",call.=FALSE)
		}
		if(chipname=="" & object@chip=="" & is.null(cdfname)){
			entrez<-refseq<-symbol<-omim<-ug<-fullname<-FALSE
			warning("Since the chip type has been specified neither",
				" by the ",type," object\n","nor by 'chipname' or",
				" 'cdfname', 'entrez', 'refseq', 'symbol',\n",
				"'omim', 'ug' and 'fullname' are set to FALSE.",call.=FALSE)
		}
	}
	if(any(c(entrez,refseq,symbol,omim,ug,fullname)))
		chipname<-check.chipname(chipname,object@chip,cdfname)
	msg<-object@msg
	h2<-unlist(strsplit(msg[1],"\n"))[1]
	h2<-substring(h2,ifelse(isSAM,13,14),nchar(h2))
	isSNP<-h2==" for Categorical Data"
	if(!is.null(refsnp)){
		tmp<-ifelse(!isSAM,"z","d")
		if(is.null(names(slot(object,tmp)))){
			refsnp<-NULL
			warning("Since no SNP names are specified by 'object'",
				"'refsnp' is ignored.",call.=FALSE)
		}
	}
	suffix<-tolower(substring(filename,nchar(filename)-4,nchar(filename)))
	if(suffix!=".html"){
		filename<-paste(filename,"html",sep=".")
		warning("Since the suffix of 'filename' is not 'html' '.html' is added",
			" to 'filename'.",call.=FALSE)
	}
	if(is.null(plotname))
		plotname<-paste(type,"plot_for_",gsub(".html","",basename(filename)),
			".png",sep="")
	file.sep<-.Platform$file.sep
	tmp<-unlist(strsplit(filename,file.sep))
	path<-if(length(tmp)>1) paste(tmp[-length(tmp)],collapse=file.sep)
		else "."	
	bg.col<-col2hex(bg.col)
	text.col<-col2hex(text.col)
	link.col<-col2hex(link.col)
	sum.out<-summary(object,delta,n.digits=n.digits,bonf=bonf,chip=chipname)
	mat.fdr<-pretty.mat.fdr(sum.out@mat.fdr,digits=n.digits)
	mat.sig<-sum.out@mat.sig
	if(!is.null(mat.sig))
		mat.sig<-pretty.mat.sig(mat.sig,digits=n.digits)
	else{
		if(addGenes){
			addGenes<-FALSE
			warning("Since there are no significant genes, 'addGenes' is set",
				" to FALSE.",call.=FALSE)
		}
	}
	if(is.null(varName))
		varName<-ifelse(isSNP,"SNPs","Genes")
	tmp.ids<-which(colnames(mat.sig)=="stdev")
	mat.sig<-mat.sig[,-tmp.ids]
	rmNA<-rowMeans(is.na(mat.sig))
	if(any(rmNA==1)){
		tmp.ids<-which(rmNA==1)
		mat.sig<-mat.sig[,-tmp.ids]
	}
	#if(isSAM && isSNP){
	#	tmp.ids<-which(colnames(mat.sig)%in%c("stdev","R.fold"))
	#	mat.sig<-mat.sig[,-tmp.ids]
	#}
	outfile<-file(filename,"w")
	cat("<html>","<head>",paste("<title>",type," Analysis</title>",sep=""),"</head>",
		paste("<body bgcolor=",bg.col," text=",text.col," link=",link.col,
		">",sep=""),paste("<h1 align=center>",type,"Analysis </h1>"),
		"<style type=text/css>",
		"h3{ margin-top: 5px; margin-bottom: 2px; padding-left: 10px; text-indent: 10px;
			word-spacing: 1px}", 
		"li{ margin-top: 10px; padding-left: 0px; text-indent: 5px; word-spacing: 1px}",
		"ul{ padding-left: 50px}",	
		"</style>",
		paste("<h2 align=center>",h2,"</h2>"),
		sep="\n",file=outfile)
	if(addPlot | !is.null(findA0)){
		suf.plot<-unlist(strsplit(plotname,"\\."))
		suf.plot<-suf.plot[length(suf.plot)]
		if(!suf.plot%in%c("jpeg","png"))
			stop("'plotname' must be either a png or a jpeg file.")
		tmp.local<-length(unlist(strsplit(plotname,file.sep)))
		plotname2<-if(tmp.local==1) paste(path,plotname,sep=file.sep) else plotname
	}
	if(!is.null(findA0)){
		pn.finda0<-gsub(type,"finda0",c(plotname,plotname2))
		finda02html(findA0,delta,outfile,plotArgs=plotFindArgs,plotnames=pn.finda0,
			bg.plot.adjust=bg.plot.adjust,bg.col=bg.col,n.digits=n.digits,
			tableborder=tableborder,plotborder=plotborder)
	}
	if(addStats){
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat(paste("<h3>General Information",
			if(!isSAM & !is.null(findA0)) " for the Actual EBAM Analysis","</h3>",
			sep=""),sep="\n",file=outfile)
		cat("<ul>",sep="\n",file=outfile)
		s0<-slot(object,ifelse(isSAM,"s0","a0"))
		if(length(s0)==1){
			s0.msg<-if(isSAM) msg[substring(msg,1,2)=="s0"] else round(s0,n.digits)
			if(length(s0.msg)==1)
				cat("<li>Fudge Factor: ",s0.msg,"\n",sep="",file=outfile)
		}
		cat(paste("<li>Prior Probability: p0 = ",
			ifelse(isSAM,mat.fdr[,"p0"],round(object@p0,n.digits)),sep=""),
			paste("<li>Statistics for Delta = ",mat.fdr[,"Delta"],":",sep=""),
			"<ul type=\"disc\">",
			paste("<li>Number of Identified ",varName,": ",
			ifelse(isSAM,mat.fdr[,"Called"],mat.fdr[,"Number"]),sep=""),
			if(isSAM)  paste("<li>Number of Falsely Called ",varName,": ",mat.fdr[,"False"],sep=""),
			paste("<li>False Discovery Rate: ",mat.fdr[,"FDR"],sep=""),
			"</ul>","</ul>",sep="\n",file=outfile)
	}
	if(addPlot){
		FUN<-match.fun(suf.plot)
		FUN(plotname2)
		if(bg.plot.adjust)
			par(bg=bg.col)
		else
			par(bg="white")
		if(isSAM)
			plot(object,delta,pos.stats=plotArgs$pos.stats,sig.col=plotArgs$sig.col,
				xlim=plotArgs$xlim,ylim=plotArgs$ylim,main=plotArgs$main,
				xlab=plotArgs$xlab,ylab=plotArgs$ylab,pty=plotArgs$pty,
				lab=plotArgs$lab,pch=plotArgs$pch,sig.cex=plotArgs$sig.cex,...)
		else
			plot(object,delta,pos.stats=plotArgs$pos.stats,sig.col=plotArgs$sig.col,
				sig.cex=plotArgs$sig.cex,pch=plotArgs$pch,stats.cex=plotArgs$stats.cex,
				main=plotArgs$main,xlab=plotArgs$xlab,ylab=plotArgs$ylab,
				y.intersp=plotArgs$y.intersp,...)
		dev.off()
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<div style=\"text-align: center\"><img src=",
			if(plotname==plotname2) "file:///",plotname," border=",
			plotborder,"></div>","\n",sep="",file=outfile)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
	}
	if(addGenes){
		has.nonames<-all(rownames(mat.sig)==as.character(1:nrow(mat.sig)))
		if(has.nonames)
			tr<-make.tablecode(as.character(mat.sig[,"Row"]),entrez=FALSE,
				refseq=FALSE,symbol=FALSE,omim=FALSE,ug=FALSE,fullname=FALSE,
				chipname=chipname,dataframe=mat.sig[,-1],load=load,
				new.window=new.window,tableborder=tableborder,name1stcol="Row")
		else
			tr<-make.tablecode(rownames(mat.sig),entrez=entrez,refseq=refseq,
				symbol=symbol,omim=omim,ug=ug,fullname=fullname,
				chipname=chipname,cdfname=cdfname,dataframe=mat.sig[,-1],
				new.window=new.window,tableborder=tableborder,refsnp=refsnp,
				which.refseq=which.refseq,load=load,max.associated=max.associated)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<h3 align=center>","Identified ",varName," (Using Delta = ",
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
		cat("The ",type," Plot required by the html file is stored in ",plotname2,
			".\n",sep="")
	if(!is.null(findA0))
		cat("The FindA0 plot required by the html file is stored in ",pn.finda0[2],
			".\n",sep="")
} 

