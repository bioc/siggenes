require(methods)
setClass("SAM",representation(d="numeric",d.bar="numeric",vec.false="numeric",p.value="numeric",
	s="numeric",s0="numeric",mat.samp="matrix",p0="numeric",mat.fdr="matrix",q.value="numeric",
	fold="numeric",msg="character",chip="character"))

setMethod("show","SAM",
	function(object){
		cat(object@msg[1])
		print(round(object@mat.fdr[,1:5],3))
	}
)

setMethod("print","SAM",
	function(x,delta=NULL,n.digits=3){
		cat(x@msg[1])
		if(is.null(delta))
			print(round(x@mat.fdr[,1:5],n.digits))
		else{
			mat.fdr<-stats.cal(x@d,x@d.bar,x@vec.false,x@p0,delta=delta)
			print(round(mat.fdr[,1:5],n.digits))
		}
	}
)


setMethod("summary","SAM",
	function(object,delta=NULL,n.digits=5,what="both",ll=TRUE,chip="",file="",sep="\t",
			quote=FALSE){
		if(is.null(delta)){
			cat(object@msg)
			print(round(object@mat.fdr,n.digits))
			return(invisible(object@mat.fdr))
		}
		if(length(delta)>1){
			cat(object@msg)
			print(round(mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
				delta=delta),n.digits))
			return(invisible(mat.fdr))
		}
		if(length(delta)==1){
			if(!what%in%c("both","stats","genes"))
				stop("what must be either stats, genes or both.")
			mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
				delta=delta)
			mat.out<-round(mat.fdr,n.digits)
			cat(object@msg[1],file=file)
			if(what%in%c("stats","both")){
				cat(object@msg[-1],file=file,append=TRUE)
				cat(" Delta:",mat.out[,"Delta"],"\n","cutlow:",mat.out[,"cutlow"],"\n",
				"cutup:",mat.out[,"cutup"],"\n","p0:",mat.out[,"p0"],"\n",
				"Significant Genes:",mat.out[,"Called"],"\n","Falsely Called Genes:",
				mat.out[,"False"],"\n","FDR:",mat.out[,"FDR"],
				ifelse(what=="both","\n\n\n","\n"),file=file,append=TRUE)
			}
			sig.genes<-which(object@d>=mat.fdr[,"cutup"] | object@d<=mat.fdr[,"cutlow"])
			mat.sig<-NULL
			if(what%in%c("genes","both")){
				if(chip=="" & object@chip=="" & ll){
					ll<-FALSE
					warning("Since the chip type is neither specified by 'chip' ",
						"nor by the SAM object,","\n","ll is set to FALSE.",
						call.=FALSE)
				}
				mat.sig<-cbind(Row=1:length(object@d),d.value=object@d,stdev=object@s,
					p.value=object@p.value,q.value=object@q.value,
					R.fold=object@fold)
				mat.sig<-mat.sig[sig.genes,]
				mat.sig<-mat.sig[rev(order(abs(mat.sig[,"d.value"]))),]
				mat.out<-as.data.frame(cbind(round(mat.sig,n.digits),
					Name=names(object@d)[mat.sig[,"Row"]]))
				mat.sig<-as.data.frame(mat.sig)
				if(ll){
					if(chip=="")
						chip<-object@chip
					if(chip!=object@chip & object@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-getLL(as.character(mat.out[,"Name"]),chip)
					sym<-getSYMBOL(as.character(mat.out[,"Name"]),chip)
					mat.sig<-data.frame(Symbol=sym,LocusLink=LL,mat.sig)
					mat.out<-data.frame(Symbol=sym,LocusLink=as.character(LL),
						mat.out)	
				}
				dimnames(mat.out)[[1]]<-1:nrow(mat.out)
				cat("Genes called significant:\n\n",file=file,append=TRUE)
				if(file=="")
					print(mat.out)
				else{
					write.table(t(dimnames(mat.out)[[2]]),file=file,sep=sep,
						append=TRUE,row.names=FALSE,col.names=FALSE,
						quote=quote)
					write.table(mat.out,file=file,sep=sep,append=TRUE,
						row.names=FALSE,col.names=FALSE,quote=quote)
				}
			}
			if(file!="")
				cat("Output is stored in",file,"\n")		 
		invisible(list(row.sig.genes=sig.genes,mat.fdr=mat.fdr,mat.sig=mat.sig))
		}
			
	}
)

setMethod("plot","SAM",
	function(x,y,pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,ylab=NULL,
		pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,helplines=TRUE,...){
		if(missing(y)){
			y<-NULL
			cat("To obtain a SAM plot, delta has to be specified in plot(object,delta,...).",
				"\n")
		}	
		if(length(y)==1)
			sam.plot2(x,delta=y,pos.stats=pos.stats,sig.col=sig.col,xlim=xlim,ylim=ylim,
				main=main,xlab=xlab,ylab=ylab,pty=pty,lab=lab,pch=pch,sig.cex=sig.cex,
				...)
		else
			delta.plot(x,delta=y,helplines=helplines)
	}
)

setMethod("identify","SAM",
	function(x,showText=TRUE,getInfo=TRUE,pos=4,cex=0.8,add.xy=numeric(2),n.digits=4,
			ask=TRUE,ll=TRUE,browse=FALSE,chip="",...){
		if(length(add.xy)!=2)
			stop("add.xy must have length 2.")
		d<-x@d
		if(is.null(names(d)))
			names(d)<-paste("G",1:length(d),sep="")
		sorted.d<-sort(d,index.return=TRUE)
		d.sort<-sorted.d$x
		mat.xy<-cbind(x@d.bar,d.sort)
		if(getInfo & chip=="" & x@chip=="" & ll){
			ll<-FALSE
			warning("Since the chip type is neither specified by 'chip' nor by the ",
				"SAM object,","\n","ll is set to FALSE.",call.=FALSE)
			}
		repeat{
			id.out<-identify(mat.xy,labels=names(d.sort),pos=TRUE,n=1,plot=FALSE)
			if(length(id.out$pos)==0)
				break
			ind<-id.out$ind
			x.ided<-mat.xy[ind,1]
			y.ided<-mat.xy[ind,2]
			if(showText)
				text(x.ided+add.xy[1],y.ided+add.xy[2],names(d.sort)[ind],pos=pos,
					cex=cex,...)
			if(any((x.ided-mat.xy[-ind,1])^2+(y.ided-mat.xy[-ind,2])^2==0)){
				same.ided<-which((x.ided-mat.xy[,1])^2+(y.ided-mat.xy[,2])^2==0)
				if(!getInfo)
					warning("There are ",length(same.ided)-1," other genes ",
						"having the same coordinates as ",names(d.sort)[ind],
						".",call.=FALSE)
			}
			else
				same.ided<-ind
			if(getInfo){
				which.d<-sorted.d$ix[same.ided]
				mat.ided<-cbind(d.value=na.exclude(d)[which.d],
					stdev=if(length(x@s)==0) numeric(0) 
						else na.exclude(x@s)[which.d],
					p.value=na.exclude(x@p.value)[which.d],
					q.value=na.exclude(x@q.value)[which.d],
					R.fold=x@fold[!is.na(d)][which.d])
				mat.ided<-data.frame(round(mat.ided,n.digits))
				if(ll){
					if(chip=="")
						chip<-x@chip
					if(chip!=x@chip & x@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-getLL(row.names(mat.ided),chip)
					sym<-getSYMBOL(row.names(mat.ided),chip)
					mat.ided<-data.frame(Symbol=sym,LocusLink=LL,mat.ided)
					if(browse){
						for(i in 1:length(na.exclude(LL))){
							browseURL(getQuery4LL(na.exclude(LL)[i]))
							if(i<length(na.exclude(LL))){
								answer2<-readline("Next LocusLink?")
								if(tolower(answer2)%in%c("n","no"))
									break
							}
						}
					}
				}		
				print(mat.ided)
				cat("\n")
				if(ask){
					answer<-readline("Press <Return> to continue. Type 'n' to stop. ")
					if(tolower(answer)%in%c("n","no"))
						break
				}
				cat("\n")				
			}
		}
        invisible(id.out)
	}
)
 
 

