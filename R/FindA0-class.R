require(methods)

setClass("FindA0",representation(mat.z="matrix",mat.posterior="matrix",mat.center="matrix",
	mat.success="matrix",mat.failure="matrix",z.norm="numeric",p0="numeric",
	mat.a0="data.frame",mat.samp="matrix",vec.a0="numeric",suggested="numeric",
	delta="numeric",df.ratio="numeric",msg="character",chip="character"))

setMethod("show","FindA0",function(object) print(object))

setMethod("print","FindA0",
	function(x,delta=NULL){
		if(is.null(delta))
			delta<-x@delta
		if(length(delta)>1)
			stop("In find.a0, the statistics can be computed only for one value of delta.")
		if(delta<=0 | delta>1)
			stop("delta must be between 0 and 1.")
		cat(x@msg[1])
		cat("Selection Criterion: Posterior >=",delta,"\n\n")
		if(delta==x@delta){
			print(x@mat.a0)
			sugg<-x@suggested
		}
		else{
			tmp<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,
				nrow(x@mat.samp),delta=delta)
			print(tmp$tab)
			sugg<-tmp$suggest
		}
		cat("\n","Suggested Choice for a0: ",round(sugg,4),sep="")
		if(names(sugg)!="-")
			cat(" (the ",100*as.numeric(names(sugg)),"% quantile of the s values)",
				sep="")
		cat("\n\n")
	}
)

setMethod("plot","FindA0",
	function(x,y,logit=TRUE,pos.legend=NULL,legend.cex=0.8,col=NULL,main=NULL,xlab=NULL,
			ylab=NULL,only.a0=FALSE,lty=1,lwd=1,y.intersp=1.1,...){
		mat.post<-x@mat.posterior
		z<-x@z.norm
		if(logit)
			mat.post<-log(mat.post/(1-mat.post))
		if(missing(y))
			y<-x@delta
		if(is.null(main))
			main<-paste("Transformed z Values vs. ",if(logit) "Logit of the ", 
				"Posterior",sep="")
		if(is.null(xlab))
			xlab<-"Transformed z Value"
		if(is.null(ylab))
			ylab<-ifelse(logit,"logit(Posterior)","Posterior")
		if(is.null(col))
			col<-1:ncol(mat.post)
		if(!length(col)%in%c(1,ncol(mat.post)))
			stop("col must be either of length 1 or equal to the number of a0 values.")
		ids<-is.finite(mat.post)
		if(any(!ids))
			warning("Some of the logit posterior probabilites are Inf. ",
				"These probabilities are not plotted.",call.=FALSE)
		ylim<-c(0,max(mat.post[ids]))
		if(is.null(pos.legend))
			pos.legend<-ifelse(any(z<0),1,4)
		if(!pos.legend%in%(0:4))
			stop("pos.legend must be an integer between 0 and 4.")
		plot(range(z),ylim,type="n",main=main,xlab=xlab,ylab=ylab)
		for(i in 1:ncol(mat.post))
			lines(z,mat.post[,i],col=col[i],lwd=lwd,lty=lty)
		h<-ifelse(logit,log(y/(1-y)),y)
		abline(h=h,lty="dashed")
		if(pos.legend!=0){
			if(y!=x@delta)
				mat.legend<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,
					nrow(x@mat.samp),delta=y)$tab
			else
				mat.legend<-x@mat.a0
			textLegend<-round(mat.legend[,1],3)
			titleLegend<-"a0"
			if(!only.a0){
				textLegend<-paste(textLegend," (",mat.legend[,3],")",sep="")
				titleLegend<-"a0 (Number)"
			}
			where<-switch(pos.legend,"top","bottomright","bottomleft","topleft")
			legend(where,legend=textLegend,lty=lty,lwd=lwd,cex=legend.cex,col=col,
				bty="n",title=titleLegend,y.intersp=y.intersp)
		}
	}
)
			




			


					