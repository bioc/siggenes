finda02html<-function(x,delta,file,plotArgs=plotFindArguments(),plotnames=NULL,
		bg.plot.adjust=FALSE,bg.col=col2hex("white"),n.digits=4,tableborder=1,
		plotborder=0){
	if(!is(x,"FindA0"))
		stop("findA0 must be an object of class FindA0.")
	if(length(plotnames)!=2)
		stop("plotnames must have length 2.")
	if(delta==x@delta){
		mat.a0<-x@mat.a0
		sugg<-x@suggested
	}
	else{
		tmp<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,nrow(x@mat.samp),
			delta=delta)
		mat.a0<-tmp$tab
		sugg<-tmp$suggest
	}
	mat.a0[,"FDR"]<-formatSAM(mat.a0[,"FDR"],digits=n.digits)
	mat.a0<-format(mat.a0,digits=n.digits)
	cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
	cat("<h3>Specification of the Fudge Factor</h3>",sep="\n",file=file)
	codeA0<-make.tablecode(rownames(mat.a0),ll=FALSE,refseq=FALSE,symbol=FALSE,omim=FALSE,
		ug=FALSE,dataframe=mat.a0,new.window=FALSE,tableborder=tableborder,
		name1stcol="&nbsp;")
	cat("<p><font color=",bg.col," size=1> HALLO</font></p>","\n",sep="",file=file)
	cat("<style type=text/css>","p{ margin-top: 2px; margin-bottom: 2px; word-spacing: 1px}",
		"</style>",codeA0,sep="\n",file=file)
	if(!plotArgs$onlyTab){
		suf.plot<-unlist(strsplit(plotnames[1],"\\."))
		suf.plot<-suf.plot[length(suf.plot)]
		FUN<-match.fun(suf.plot)
		FUN(plotnames[2])
		if(bg.plot.adjust)
			par(bg=bg.col)
		else
			par(bg="white")
		plot(x,delta,logit=plotArgs$logit,pos.legend=plotArgs$pos.legend,
			cexLegend=plotArgs$cexLegend,col=plotArgs$col,main=plotArgs$main,
			xlab=plotArgs$xlab,ylab=plotArgs$ylab,only.a0=plotArgs$only.a0,
			lty=plotArgs$lty,lwd=plotArgs$lwd,y.intersp=plotArgs$y.intersp)
		dev.off()
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
		cat("<div style=\"text-align: center\"><img src=",
			if(plotnames[1]==plotnames[2]) "file:///",plotnames[1]," border=",
			plotborder,"></div>","\n",sep="",file=file)
		cat("<p><b>Suggested Choice:</b>"," a0 = ",round(sugg,n.digits),
			"   (Selection Criterion: Number of Genes with Posterior >= ",delta,")",
			"</p>","\n",sep="",file=file)
	}
	cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
}
		
		

	