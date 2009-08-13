denspr<-function(x,n.interval=NULL,df=5,knots.mode=TRUE,type.nclass=c("wand","scott","FD"),
		addx=FALSE){
	require(splines)
	if(is.null(n.interval)){
		type<-match.arg(type.nclass)
		FUN<-match.fun(paste("nclass",type,sep="."))
		n.interval<-FUN(x)
	}
	else
		type<-NULL
	breaks<-seq(min(x),max(x),length=n.interval+1)
	valHist<-hist(x,breaks=breaks,plot=FALSE)
	center<-valHist$mids
	counts<-valHist$counts
	ids<-which(counts>0)
	x.mode<-center[which.max(counts)]
	if(knots.mode){
		x.q<-mean(center<=x.mode)
		q.knots<-getQuantiles(df-1,x.q)
	}
	center<-center[ids]
	if(knots.mode){
		knots<-quantile(center,q.knots)
		tmp<-ns.out<-ns(center,knots=knots)
	}
	else
		tmp<-ns.out<-ns(center,df=df)
	class(tmp)<-"matrix"
	mat<-data.frame(Number=counts[ids],tmp)
	glm.out<-glm(Number~.,data=mat,family=poisson)
	scale<-sum(diff(breaks)*counts)
	newx<-predict(ns.out,x)
	class(newx)<-"matrix"
	preds<-predict(glm.out,data.frame(newx),type="response")
	out <- list(y=preds/scale, center=valHist$mids, counts=counts, x.mode=x.mode,
		ns.out=ns.out, type=type, x=if(addx) x else NULL)
	class(out) <- "denspr"
	out
}


plot.denspr <- function(x, ylab="Density", xlab="x", type="l", ...){
	if(is.null(x$x))
		stop("If the density should be plotted, addx must be set to TRUE in denspr.")
	xval <- unique(x$x)
	y <- unique(x$y)
	plot(xval[order(xval)], y[order(xval)], type=type, ylab=ylab, xlab=xlab, ...)
}


lines.denspr <- function(x, type="l", ...){
	if(is.null(x$x))
		stop("If the density should be plotted, addx must be set to TRUE in denspr.")
	xval <- unique(x$x)
	y <- unique(x$y)
	lines(xval[order(xval)], y[order(xval)], type=type, ...)
}

