plotArguments<-function(pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,
		ylab=NULL,pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,cexStats=0.8,
		y.intersp=1.3){
	list(pos.stats=pos.stats,sig.col=sig.col,xlim=xlim,ylim=ylim,main=main,
		xlab=xlab,ylab=ylab,pty=pty,lab=lab,pch=pch,sig.cex=sig.cex,
		cexStats=cexStats,y.intersp=y.intersp)
}


plotFindArguments<-function(onlyTab=FALSE,logit=TRUE,pos.legend=NULL,cexLegend=0.8,col=NULL,
		main=NULL,xlab=NULL,ylab=NULL,only.a0=FALSE,lty=1,lwd=1,y.intersp=1.1){
	list(onlyTab=onlyTab,logit=logit,pos.legend=pos.legend,cexLegend=cexLegend,col=col,
		main=main,xlab=xlab,ylab=ylab,only.a0=only.a0,lty=lty,lwd=lwd,
		y.intersp=y.intersp)
}
