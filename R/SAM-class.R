require(methods)
setClass("SAM",representation(d="numeric",d.bar="numeric",vec.false="numeric",p.value="numeric",
	s="numeric",s0="numeric",mat.samp="matrix",p0="numeric",mat.fdr="matrix",q.value="numeric",
	fold="numeric",msg="character",chip="character"))

setMethod("show","SAM",
	function(object){
		cat(object@msg[1])
		print(pretty.mat.fdr(object@mat.fdr[,1:5],digits=3))
	}
)

setMethod("print","SAM",
	function(x,delta=NULL,n.digits=3){
		cat(x@msg[1])
		if(is.null(delta))
			print(pretty.mat.fdr(x@mat.fdr[,1:5],digits=n.digits))
		else{
			mat.fdr<-as.data.frame(stats.cal(x@d,x@d.bar,x@vec.false,x@p0,
				delta=delta))
			print(pretty.mat.fdr(mat.fdr[,1:5],digits=n.digits))
		}
	}
)



setMethod("summary","SAM",
	function(object,delta=NULL,n.digits=3,what="both",entrez=FALSE,bonf=FALSE,
			chip="",file="",sep="\t",quote=FALSE,dec="."){
		list.args<-list(n.digits=n.digits,what=what,file=file,sep=sep,quote=quote,
			dec=dec,msg=object@msg)
		if(length(delta)!=1){
			mat.fdr<-if(is.null(delta)) object@mat.fdr
				else  stats.cal(object@d,object@d.bar,object@vec.false,
					object@p0,delta=delta)
			sig.genes<-numeric(0)
			mat.sig<-data.frame(NULL)
		}
		else{
			if(!what%in%c("both","stats","genes"))
				stop("'what' must be either \"stats\", \"genes\" or \"both\".")
			mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
				delta=delta)
			sig.genes<-which(object@d>=mat.fdr[,"cutup"] | object@d<=mat.fdr[,"cutlow"])
			mat.sig<-data.frame(NULL)
			if(what%in%c("genes","both") & length(sig.genes)!=0){
				if(chip=="" & object@chip=="" & entrez){
					entrez<-FALSE
					warning("Since the chip type is neither specified by 'chip' ",
						"nor by the SAM object,","\n","entrez is set to FALSE.",
						call.=FALSE)
				}
				if(bonf)
					p.bonf <- p.adjust(object@p.value, "bonferroni")
				mat.sig<-cbind(Row=(1:length(object@d))[sig.genes],
					d.value=object@d[sig.genes],stdev=object@s[sig.genes],
					rawp=object@p.value[sig.genes],
					q.value=object@q.value[sig.genes],
					Bonferroni = if(bonf) p.bonf[sig.genes], 
					R.fold=object@fold[sig.genes])
				if(length(sig.genes)>1)
					mat.sig<-mat.sig[rev(order(abs(mat.sig[,"d.value"]))),]
				mat.sig<-as.data.frame(mat.sig)
				if(entrez & all(row.names(mat.sig)==as.character(1:nrow(mat.sig)))){
					entrez<-FALSE
					warning("Since no gene names are available it is not",
						" possible to obtain Entrez links.\n",
						"'entrez' is thus set to FALSE.",call.=FALSE)
				}
				if(entrez){
					if(chip=="")
						chip<-object@chip
					if(chip!=object@chip & object@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-unlist(lookUp(row.names(mat.sig),chip,"ENTREZID"))
					sym<-getSYMBOL(row.names(mat.sig),chip)
					mat.sig<-data.frame(Row=mat.sig[,1],Symbol=sym,
						Entrez=LL,mat.sig[,-1])
				}
			}
			retval<-new("sumSAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,
				mat.sig=mat.sig,list.args=list.args)
		}
		new("sumSAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,mat.sig=mat.sig,
			list.args=list.args)
	}
)

setMethod("plot","SAM",
	function(x,y,pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,ylab=NULL,
		pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,helplines=FALSE,...){
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
	function(x,showText=TRUE,getInfo=TRUE,pos=4,cex=0.8,add.xy=numeric(2),n.digits=3,
			ask=FALSE,entrez=FALSE,browse=FALSE,chip="",...){
		if(length(add.xy)!=2)
			stop("add.xy must have length 2.")
		d<-x@d
		if(is.null(names(d)))
			names(d)<-paste("G",1:length(d),sep="")
		sorted.d<-sort(d,index.return=TRUE)
		d.sort<-sorted.d$x
		mat.xy<-cbind(x@d.bar,d.sort)
		if(getInfo & chip=="" & x@chip=="" & entrez){
			entrez<-FALSE
			warning("Since the chip type is neither specified by 'chip' nor by the ",
				"SAM object,","\n","entrez is set to FALSE.",call.=FALSE)
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
					rawp=na.exclude(x@p.value)[which.d],
					q.value=na.exclude(x@q.value)[which.d],
					R.fold=x@fold[!is.na(d)][which.d])
				mat.ided<-data.frame(round(mat.ided,n.digits))
				if(entrez){
					if(chip=="")
						chip<-x@chip
					if(chip!=x@chip & x@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-unlist(lookUp(row.names(mat.ided),chip,"ENTREZID"))
					sym<-getSYMBOL(row.names(mat.ided),chip)
					mat.ided<-data.frame(Symbol=sym,Entrez=LL,mat.ided)
					if(browse){
						for(i in 1:length(na.exclude(LL))){
							browseURL(annotate:::getQuery4LL(na.exclude(LL)[i]))
							if(i<length(na.exclude(LL))){
								answer2<-readline("Next Entrez link?")
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



