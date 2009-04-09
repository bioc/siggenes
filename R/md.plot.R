md.plot<-function(object, delta, pos.stats=1, sig.col=3, xlim=NULL, ylim=NULL, main=NULL,
		xlab=NULL, ylab=NULL, xsym=NULL, ysym=NULL, forceDelta=FALSE, includeZero=TRUE,
		lab=c(10,10,7), pch=NULL, sig.cex=1, ...){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	#if(is.null(pos.stats))
	#	pos.stats<-ifelse(all(object@d>=0,na.rm=TRUE),2,1)
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
	if(any(d.sort<0))
		xsym <- ysym <- TRUE
	else
		xsym <- ysym <- FALSE
	M <- d.sort - object@d.bar
	if(is.null(xlim)){
		if(xsym){
			xmax <- max(abs(d.sort))
			xlim <- c(-xmax, xmax)
		}
		else
			xlim<-range(d.sort, if(includeZero) 0)
	}
	if(is.null(ylim)){
		if(ysym){
			ymax <- max(abs(M))
			ylim <- c(-ymax, ymax)
		}
		else
			ylim<-range(M, if(includeZero) 0)
		if(forceDelta)
			ylim <- range(ylim, delta, -delta)
	}
	if(is.null(main))
		main<-paste("MD Plot for Delta =",delta)
	if(is.null(xlab))
		xlab<-"D"
	if(is.null(ylab))
		ylab<-"M"
	par.store<-par()$lab
	on.exit(par(lab=par.store))
	par(lab=lab)
	mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,delta=delta)
	d.up<-which(d.sort>=mat.fdr[,"cutup"])
	d.down<-which(d.sort<=mat.fdr[,"cutlow"])
	stats<-paste(c("cutlow:","cutup:","p0:","Significant:","False:","FDR:"),
		round(mat.fdr[1,c("cutlow","cutup","p0","Called","False","FDR")],3),sep="  ")
	plot(xlim, ylim, type="n", xlab=xlab, ylab=ylab, main=main)
	plotLegendLines(stats, pos.stats, delta, mat.fdr[,"cutup"], mat.fdr[,"cutlow"], par("usr"))
	if(length(c(d.up,d.down))==0)
		points(d.sort,M,main=main,pch=pch,...)
	else{
		points(d.sort[-c(d.up,d.down)],M[-c(d.up,d.down)],pch=pch,...)
		if(length(d.up)>0)
			points(d.sort[d.up],M[d.up],col=col.up,cex=sig.cex,pch=pch)
		if(length(d.down)>0)
			points(d.sort[d.down],M[d.down],col=col.down,cex=sig.cex,pch=pch)
	}
}

 
plotLegendLines <- function(stats, pos.stats, delta, up, low, lims){
	whereStats <- ifelse(pos.stats==1, "topleft", "bottomright")
	xlim2 <- lims[1:2]
	ylim2 <- lims[3:4]
	if(pos.stats<2){
		abline(h=-delta, lty=2)
		abline(v=up, lty=5, cex=1.5)
	}
	if(pos.stats!=1){
		abline(h=delta, lty=2)
		abline(v=low, lty=5, cex=1.5)
	}
	if(pos.stats>0){
		rect <- legend(whereStats, legend=stats, bty="n", cex=0.75, y.intersp=1.1)$rect
		top <- rect$top
		left <- rect$left
		right <- left + rect$w
		bottom <- top - rect$h
	}
	if(pos.stats==1){
		segments(ifelse(bottom<delta, right, xlim2[1]), delta, xlim2[2], delta, lty=2)
		segments(low, ylim2[1], low, ifelse(right>low, bottom, ylim2[2]), lty=5, cex=1.5)
	}
	if(pos.stats==2){
		segments(xlim2[1], -delta, ifelse(top>-delta, left, xlim2[2]), -delta, lty=2)
		segments(up, ifelse(up>left, top, ylim2[1]), up, ylim2[2], lty=5, cex=1.5)
	}
}

