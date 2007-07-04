`findFDR` <-
function(object,genes,delta=NULL,isSAM=TRUE,prec=6,verbose=FALSE,first=TRUE){
	if(isSAM)
		mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
			delta=delta)[,c(1,4,5)]
	else
		mat.fdr<-compNumber(object@z,object@posterior,object@p0,nrow(object@mat.samp),
			delta=delta,vec.pos=object@vec.pos,vec.neg=object@vec.neg)[,1:3]
	vec.num<-mat.fdr[,2]
	if(first){
		if(all(vec.num>genes)){
			cat("For Delta <= ",delta[length(delta)],", the number of variables is larger than ",genes,
				".\n",sep="")
			return(mat.fdr[1,])
		}
		if(all(vec.num<genes)){
			cat("For Delta >= ",delta[1],", the number of variables is smaller than ",
				genes,".\n",sep="")
			return()
		}
		if(verbose)
			cat(" Starting search between:\n","Delta",delta[1],"    Number:",vec.num[1],"\n",
				"Delta:",delta[10],"    Number:",vec.num[10],"\n")
	}
	if(any(vec.num==genes)){
		thres<-min(which(vec.num==genes))
		if(verbose)
			cat("\n\n","The threshold Delta is given by\n",sep="")
		print(round(mat.fdr[thres,,drop=FALSE],prec))
		return(mat.fdr[thres,])
	}
	thres<-max(which(vec.num>genes))
	delta.new<-round(seq(delta[thres],delta[thres+1],le=10),prec)
	if(all(delta==delta.new)){
		if(verbose)
			cat("\n\n")
		cat("The threshold seems to be at","\n")
		print(round(mat.fdr[c(thres,thres+1),],prec))
		return(mat.fdr[c(thres,thres+1),])
	}
	delta<-delta.new
	if(verbose)
		cat("\n","Now searching between:","\n","Delta",delta[1],"    Number:",vec.num[thres],"\n",
			"Delta:",delta[10],"    Number:",vec.num[thres+1],"\n")
	findFDR(object,genes,delta=delta,isSAM=isSAM,prec=prec,verbose=verbose,first=FALSE)
}

