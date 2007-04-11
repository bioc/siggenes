require(methods)

setClass("EBAM",representation(z="numeric",posterior="numeric",p0="numeric",local="numeric",
	mat.fdr="matrix",a0="numeric",mat.samp="matrix",vec.pos="numeric",vec.neg="numeric",
	msg="character",chip="character"))

setMethod("show","EBAM",function(object) print(object))

setMethod("print","EBAM",
	function(x,delta=NULL,n.digits=4){
		if(is.null(delta))
			mat.fdr<-x@mat.fdr[,1:3,drop=FALSE]
		else
			mat.fdr<-compNumber(x@z,x@posterior,x@p0,nrow(x@mat.samp),delta=delta,
				vec.pos=x@vec.pos,vec.neg=x@vec.neg)[,1:3,drop=FALSE]
		cat(x@msg[1])
		if(length(x@a0)>0)
			cat("Fudge Factor:  a0 =",round(x@a0,n.digits),"\n\n")
		print(round(mat.fdr,n.digits))
	}
)

setMethod("summary","EBAM",
	function(object,delta=NULL,n.digits=4,what="both",ll=FALSE,chip="",file="",
			sep="\t",quote=FALSE,dec="."){
		if(is.null(delta))
			stop("delta must be specified.")
		if(!what%in%c("both","stats","genes"))
			stop("'what' must be either \"stats\", \"genes\" ",
				"or \"both\".")
		mat.fdr<-compNumber(object@z,object@posterior,object@p0,nrow(object@mat.samp),
			delta=delta,vec.pos=object@vec.pos,vec.neg=object@vec.neg)
		sig.genes<-which(object@z>=mat.fdr[,"CU"] | object@z<=mat.fdr[,"CL"])
		if(what%in%c("genes","both") & length(sig.genes)!=0){
			mat.sig<-cbind(Row=sig.genes,z.value=object@z[sig.genes],
				posterior=object@posterior[sig.genes],
				local.fdr=object@local[sig.genes])
			row.names(mat.sig)<-names(object@z)[sig.genes]
			mat.sig<-mat.sig[rev(order(abs(mat.sig[,"z.value"]))),,drop=FALSE]
			mat.sig<-as.data.frame(mat.sig)
			if(ll){
				if(chip=="" & object@chip==""){
					ll<-FALSE
					warning("Since the chip type is neither specified by ",
						"'chip' nor by the EBAM object,\n",
						"ll is set to FALSE.",call.=FALSE)
				}
				if(all(row.names(mat.sig)==as.character(1:nrow(mat.sig)))){
					ll<-FALSE
					warning("Since no gene names are available, it is not",
						" possible to obtain locus links.\n",
						"Thus, 'll' is set to FALSE.",call.=FALSE)
				}
			}
			if(ll){
				if(chip=="")
					chip<-object@chip
				if(chip!=object@chip & object@chip!="")
					stop("'chip' differs from the chip type of the EBAM object.")
				require(annotate)
				LL<-getLL(row.names(mat.sig),chip)
				sym<-getSYMBOL(row.names(mat.sig),chip)
				mat.sig<-data.frame(Row=mat.sig[,1],Symbol=sym,LocusLink=LL,
					mat.sig[,-1])
			} 
		}
		else
			mat.sig<-data.frame(NULL)
		list.args<-list(n.digits=n.digits,what=what,file=file,sep=sep,quote=quote,
			dec=dec,msg=object@msg,p0=object@p0,a0=object@a0)
		new("sumEBAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,mat.sig=mat.sig,
			list.args=list.args)
	}
)
		
setMethod("plot","EBAM",
	function(x,y,pos.stats=NULL,sig.col=3,sig.cex=1,pch=NULL,stats.cex=0.8,main=NULL,
			xlab=NULL,ylab=NULL,y.intersp=1.3,...){
		z<-x@z
		post<-x@posterior
		if(missing(y))
			stop("No delta value has been specified.")
		if(length(y)!=1)
			stop("More than one delta value has been specified.")
		mat.fdr<-compNumber(x@z,x@posterior,x@p0,nrow(x@mat.samp),delta=y,
			vec.pos=x@vec.pos,vec.neg=x@vec.neg)
		if(is.null(main))
			main<-paste("EBAM Plot for Delta =",y)
		if(is.null(xlab))
			xlab<-"z Value"
		if(is.null(ylab))
			ylab<-"Posterior"
		if(length(sig.col)>1)
			stop("sig.col must be of length 1.")
		ids<-which(z<=mat.fdr[,4] | z>=mat.fdr[,5])
		twosided<-any(z<0)
		if(is.null(pos.stats))
			pos.stats<-ifelse(twosided,2,4)
		if(!pos.stats%in%(0:4))
			stop("pos.stats must be an integer between 0 and 4.")
		if(length(ids)==0)
			plot(z,post,main=main,xlab=xlab,ylab=ylab,pch=pch,...)
		else{
			plot(z[-ids],post[-ids],main=main,xlab=xlab,ylab=ylab,pch=pch,
				xlim=range(z),ylim=range(post),...)
			points(z[ids],post[ids],cex=sig.cex,col=sig.col,pch=pch)
		}
		abline(h=y,lty="dashed")
		if(pos.stats!=0){
			tmp<-c("Significant:","FDR:","p0:",if(length(x@a0)==1) "a0:",
				if(twosided) "Cutlow:","Cutup:")
			tmp2<-c(mat.fdr[,2],round(mat.fdr[,3],3),round(x@p0,3),
				if(length(x@a0)==1) round(x@a0,3), 
				if(twosided) round(mat.fdr[,4],3), round(mat.fdr[,5],3))
			textLegend<-paste(tmp,tmp2,sep="  ")
			where<-switch(pos.stats,"top","bottomright","bottomleft","topleft")
			legend(where,legend=textLegend,cex=stats.cex,bty="n",y.intersp=y.intersp)
		}
	}
)



			


					