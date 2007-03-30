require(methods)

setClass("sumEBAM",representation(row.sig.genes="numeric",mat.fdr="matrix",
	mat.sig="data.frame",list.args="list"))

setMethod("show","sumEBAM",function(object) print(object))

setMethod("print","sumEBAM",
	function(x){
		list.args<-x@list.args
		rsg<-x@row.sig.genes
		mat.fdr<-x@mat.fdr
		mat.sig<-x@mat.sig
		nd<-list.args$n.digits
		file<-list.args$file
		cat(list.args$msg,"Delta: ",mat.fdr[,"Delta"],"\n",sep="",file=file)
		if(length(list.args$a0)!=0)
			cat("a0: ",round(list.args$a0,nd),"\n",sep="",file=file,append=TRUE)
		cat("p0: ",round(list.args$p0,nd),"\n","Cutlow: ",round(mat.fdr[,"CL"],nd),"\n",
			"Cutup: ",round(mat.fdr[,"CU"],nd),"\n","Identified Genes: ",
			round(mat.fdr[,"Number"],nd),"\n","Estimated FDR: ",
			formatSAM(mat.fdr[,"FDR"],nd),"\n",file=file,append=TRUE,sep="")
		if(list.args$what%in%c("genes","both") & length(rsg)!=0){
			if(all(row.names(mat.sig)!=as.character(1:nrow(mat.sig))))
				mat.sig<-data.frame(mat.sig,Name=rownames(mat.sig))
			mat.sig[,"local.fdr"]<-formatSAM(mat.sig[,"local.fdr"],digits=nd)
			mat.sig<-format(mat.sig,digits=nd)
			row.names(mat.sig)<-1:nrow(mat.sig)
			cat("\n\n","Genes called differentially expressed (posterior >= ",
				mat.fdr[,"Delta"],"):\n\n",file=file,append=TRUE,sep="")
			if(file=="")
				print(mat.sig)
			else{
				write.table(t(dimnames(mat.sig)[[2]]),file=file,
					sep=list.args$sep,append=TRUE,row.names=FALSE,
					col.names=FALSE,quote=list.args$quote,dec=list.args$dec)
				write.table(mat.sig,file=file,sep=list.args$sep,append=TRUE,
					row.names=FALSE,col.names=FALSE,quote=list.args$quote,
					dec=list.args$dec)
			}
		}
		if(file!="")
			cat("Output is stored in", file,"\n")
	}
)


			


					