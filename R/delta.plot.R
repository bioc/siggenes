delta.plot<-function(object,delta=NULL,helplines=FALSE){
	par(mfrow=c(1,2))
	mat.fdr<-if(is.null(delta)) object@mat.fdr 
		else stats.cal(object@d,object@d.bar,object@vec.false,object@p0,delta=delta)
	plot(mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],main="Delta vs. FDR",xlab=expression(Delta),
		ylab="FDR (in %)",type="b")
	if(helplines){
        	segments(0,100*mat.fdr[,"FDR"],mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],lty=2)
        	segments(mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],mat.fdr[,"Delta"],
			-100*max(mat.fdr[,"FDR"]),lty=2)
    	}
	plot(mat.fdr[,"Delta"],mat.fdr[,"Called"],main="Delta vs. Significant Genes", 
        	xlab=expression(Delta),ylab="Number of Significant Genes",type="b")
    	if(helplines){
		segments(0,mat.fdr[,"Called"],mat.fdr[,"Delta"],mat.fdr[,"Called"],lty=2)
        	segments(mat.fdr[,"Delta"],mat.fdr[,"Called"],mat.fdr[,"Delta"],
			-max(mat.fdr[,"Called"]),lty=2)
	}
	par(mfrow=c(1,1))
}
 
 

