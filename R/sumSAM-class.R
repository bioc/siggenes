require(methods)

setClass("sumSAM",representation(row.sig.genes="numeric",mat.fdr="matrix",
	mat.sig="data.frame",list.args="list"))

setMethod("show","sumSAM",
	function(object) print(object)
)

setMethod("print","sumSAM",
	function(x,varNames=NULL){
		list.args<-x@list.args
		file<-list.args$file
		sep<-list.args$sep
		n.digits<-list.args$n.digits
		what<-list.args$what
		quote<-list.args$quote
		dec<-list.args$dec
		msg<-list.args$msg
		row.sig.genes<-x@row.sig.genes
		mat.fdr<-x@mat.fdr
		mat.sig<-x@mat.sig
		if(length(row.sig.genes)==0){
			cat(msg,file=file)
			mat.out<-pretty.mat.fdr(mat.fdr,digits=n.digits)
			if(file=="")
				print(mat.out)
			else{
				write.table(t(dimnames(mat.out)[[2]]),file=file,sep=sep,
					append=TRUE,row.names=FALSE,col.names=FALSE,
					quote=quote,dec=dec)
				write.table(mat.out,file=file,sep=sep,append=TRUE,
					row.names=FALSE,col.names=FALSE,quote=quote,
					dec=dec)
			}
		}
		else{	
			if(is.null(varNames)){
				h2<-unlist(strsplit(msg[1],"\n"))[1]
				h2<-substring(h2,13,nchar(h2))
				varNames<-ifelse(h2==" for Categorical Data","SNPs","Genes")
			}
			cat(msg[1],file=file)
			if(what%in%c("stats","both")){
				mat.out<-pretty.mat.fdr(mat.fdr,digits=n.digits)
				cat(msg[-1],file=file,append=TRUE)
				cat("Delta: ",mat.out[,"Delta"],"\n","cutlow: ",
					mat.out[,"cutlow"],"\n","cutup: ",mat.out[,"cutup"],
					"\n","p0: ",mat.out[,"p0"],"\n","Identified ",varNames,
					": ",
					mat.out[,"Called"],"\n","Falsely Called ",varNames,": ",
					mat.out[,"False"],"\n","FDR: ",mat.out[,"FDR"],
					ifelse(what=="both","\n\n\n","\n"),file=file,sep="",
					append=TRUE)
			}
			if(what%in%c("genes","both") & length(row.sig.genes)!=0){
				if(all(rownames(mat.sig)!=as.character(1:nrow(mat.sig))))
					mat.sig<-data.frame(mat.sig,Name=rownames(mat.sig))
				mat.out<-pretty.mat.sig(mat.sig,digits=n.digits)
				rownames(mat.out)<-1:nrow(mat.out)
				cat("Identified ",varNames," (using Delta = ",
					mat.fdr[,"Delta"],"):\n\n",file=file,append=TRUE,
					sep="")
				if(file=="")
					print(mat.out)
				else{
					write.table(t(dimnames(mat.out)[[2]]),file=file,
						sep=sep,append=TRUE,row.names=FALSE,
						col.names=FALSE,quote=quote,dec=dec)
					write.table(mat.out,file=file,sep=sep,append=TRUE,
						row.names=FALSE,col.names=FALSE,
						quote=quote,dec=dec)
				}
			}
		}
		if(file!="")
			cat("Output is stored in",file,"\n")
	}
)



			