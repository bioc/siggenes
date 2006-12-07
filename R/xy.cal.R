xy.cal<-function(cl,wilc=FALSE,emp=FALSE){
	ebs<-ifelse(emp,"EB","S")
	ana.type<-paste(ebs,ifelse(wilc,"AM-Wilc","AM"),sep="")
	lev<-unique(cl)
	uni.cl<-length(lev)
	uni.cl.abs<-length(unique(abs(cl)))
	if(uni.cl>2 & uni.cl!=2*uni.cl.abs)
		stop("There is something wrong with the classlabels.")
	if(uni.cl==1){
		if(wilc)
			stop("SAM-Wilc and EBAM-Wilc are not available yet for one class data.")
		cat(ana.type,"Analysis for the one-class case.","\n","\n")
		paired<-TRUE
		if(lev!=1)
			cat("Warning: Expected classlabel is 1. cl will thus be set to 1.","\n","\n")
		x<-1:length(cl)
		y<-NULL
	}
	if(uni.cl==2){
		cat(ana.type,"Analysis for the two class unpaired case.","\n","\n")
		paired<-FALSE
		if(min(lev)!=0 | max(lev)!=1){
			cat("Warning: Expected classlabels are 0 and 1. cl will thus be set to 0 and 1.","\n","\n")
			cl[which(cl==min(lev))]<-0
			cl[which(cl==max(lev))]<-1
		}
		x<-which(cl==1)
		y<-which(cl==0)
	}
	if(uni.cl==2*uni.cl.abs){
		cat(ana.type,"Analysis for the two class paired case.","\n","\n")
		paired<-TRUE
		sort.cl<-sort(cl,index=TRUE)
		if(!all(sort.cl$x==c(-uni.cl.abs:-1,1:uni.cl.abs)))
			stop("There is something wrong with the classlabels.")
		x<-sort.cl$ix[(uni.cl.abs+1):uni.cl]
		y<-sort.cl$ix[uni.cl.abs:1]
	}
	structure(list(x=x,y=y,paired=paired))
}
 
 

