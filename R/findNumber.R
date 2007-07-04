`findNumber` <-
function(object,fdr,delta=NULL,isSAM=TRUE,prec=6,verbose=FALSE,first=TRUE){
	if(isSAM)
		mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
			delta=delta)[,c(1,4,5)]
	else
		mat.fdr<-compNumber(object@z,object@posterior,object@p0,nrow(object@mat.samp),
			delta=delta,vec.pos=object@vec.pos,vec.neg=object@vec.neg)[,1:3]
	vec.fdr<-mat.fdr[,"FDR"]
	if(first){
		if(all(vec.fdr>fdr)){
			cat("For Delta <= ",delta[length(delta)],", all FDRs are larger than ",fdr,".\n",sep="")
			return()
		}
		if(all(vec.fdr<fdr)){
			cat("For Delta >= ",delta[1],", all FDRs are smaller than ",fdr,".\n",sep="")
			return(mat.fdr[1,])
		}
		if(verbose)
			cat(" Starting search between:","\n","Delta",delta[1],"    FDR:",vec.fdr[1],"\n",
				"Delta:",delta[10],"    FDR:",vec.fdr[10],"\n")
	}
	if(any(vec.fdr==fdr)){
		thres<-min(which(vec.fdr==fdr))
		if(verbose)
			cat("\n\n","The threshold Delta is given by\n",sep="")
		print(round(mat.fdr[thres,,drop=FALSE],prec))
		return(mat.fdr[thres,])
	}
	thres<-max(which(vec.fdr>fdr))
	delta.new<-round(seq(delta[thres],delta[thres+1],le=10),prec)
	if(all(delta==delta.new)){
		if(verbose)
			cat("\n\n")
		cat("The threshold seems to be at","\n")
		print(round(mat.fdr[c(thres,thres+1),],6))
		return(mat.fdr[c(thres,thres+1),])
	}
	delta<-delta.new
	if(verbose)
		cat("\n","Now searching between:","\n","Delta",delta[1],"    FDR:",vec.fdr[thres],"\n",
			"Delta:",delta[10],"    FDR:",vec.fdr[thres+1],"\n")
	findNumber(object,fdr,delta=delta,isSAM=isSAM,prec=prec,verbose=verbose,first=FALSE)
}

