qvalue.cal<-function(p,p0,version=1){
	p.na<-NULL
	p.names<-names(p)
	if(any(is.na(p))){
		p.na<-which(is.na(p))
		n.p<-length(p)
		p<-as.numeric(na.exclude(p))
	}
	if(!version%in%c(1,2))
		stop("version must be either 1 (for using the FDR) or 2 (pFDR).")
	m<-length(p)
	p.sort<-sort(p)
	qvalue<-p0*m*p.sort/(1:m)
	if(version==2)
		qvalue<-qvalue/(1-(1-ifelse(p.sort!=0,p.sort,min(p.sort[p.sort!=0],10^-15)))^m)
	qvalue[m]<-min(qvalue[m],1)
	for(i in (m-1):1)
		qvalue[i]<-min(qvalue[i],qvalue[i+1])
	qvalue[order(p)]<-qvalue
	if(!is.null(p.na)){
		qnew<-rep(NA,n.p)
		qnew[-p.na]<-qvalue
		qvalue<-qnew
	}
	names(qvalue)<-p.names
	qvalue
}
	
	 
 

