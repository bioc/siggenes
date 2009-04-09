sam.plot2<-function(object,delta,pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,
		ylab=NULL,pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,...){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	if(is.null(pos.stats))
		pos.stats<-ifelse(all(object@d>=0,na.rm=TRUE),2,1)
	if(!pos.stats%in%0:2)
		stop("pos.stats must be either 0 (statistics are not displayed),\n",
			"1 (stats are displayed in the upper left of the plot), or 2 (lower right).")
	if(length(sig.col)==1)
		col.down<-col.up<-sig.col
	else{
		col.down<-sig.col[1]
		col.up<-sig.col[2]
	}
	d.sort<-sort(object@d)
	d.bar<-object@d.bar
	if(is.null(xlim))
		xlim<-c(min(d.sort,d.bar),max(d.sort,d.bar))
	if(is.null(ylim))
		ylim<-c(min(d.sort,d.bar),max(d.sort,d.bar))
	if(is.null(main))
		main<-paste("SAM Plot for Delta =",delta)
	if(is.null(xlab))
		xlab<-"Expected d(i) values"
	if(is.null(ylab))
		ylab<-"Observed d(i) values"
	par.store<-list(par()$pty,par()$lab)
	on.exit(par(pty=par.store[[1]],lab=par.store[[2]]))
	par(pty=pty,lab=lab)
	mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,delta=delta)
	d.up<-which(d.sort>=mat.fdr[,"cutup"])
	d.down<-which(d.sort<=mat.fdr[,"cutlow"])
	if(length(c(d.up,d.down))==0)
		plot(d.bar,d.sort,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,pch=pch,...)
	else{
		plot(d.bar[-c(d.up,d.down)],d.sort[-c(d.up,d.down)],main=main,xlab=xlab,ylab=ylab,
			xlim=xlim,ylim=ylim,pch=pch,...)
		points(d.bar[d.up],d.sort[d.up],col=col.up,cex=sig.cex,pch=pch)
		points(d.bar[d.down],d.sort[d.down],col=col.down,cex=sig.cex,pch=pch)
	}
	abline(0,1)
	abline(delta,1,lty=2)
	abline(-delta,1,lty=2)
	abline(h=mat.fdr[,"cutup"],lty=5,cex=1.5)
	abline(h=mat.fdr[,"cutlow"],lty=5,cex=1.5)
	stats<-paste(c("cutlow:","cutup:","p0:","Significant:","False:","FDR:"),
		round(mat.fdr[1,c("cutlow","cutup","p0","Called","False","FDR")],3),sep="  ")
	if(pos.stats==1)
		text(rep(xlim[1],6),seq(ylim[2],ylim[2]-(ylim[2]-ylim[1])/4,le=6),stats,
			adj=0,cex=.75)
	if(pos.stats==2)
		text(rep(xlim[2]-(xlim[2]-xlim[1])/4,6),
			seq(ylim[1],ylim[1]+(ylim[2]-ylim[1])/4,le=6),stats,adj=0,cex=.75)
}

 
 

