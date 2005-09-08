stats.cal<-function(d,d.bar,vec.false,p0,delta=NULL,le.delta=10){
	d.sort<-sort(d)
	d.diff<-d.sort-d.bar
	m<-length(d.diff)
	if(is.null(delta)){
		ra.ddiff<-range(abs(d.diff))
		delta<-round(seq(max(0.1,ra.ddiff[1]),max(1,ra.ddiff[2]),le=le.delta),1)
	}
	else{
		if(any(delta<=0))
			stop("delta must be larger than 0.")
		le.delta<-length(delta)
	}
	j0<-which(d.bar==min(d.bar[d.bar>=0]))[1]
	mat.fdr<-matrix(0,le.delta,9)
	dimnames(mat.fdr)<-list(1:le.delta,
		c("Delta","p0","False","Called","FDR","cutlow","cutup","j2","j1"))
	mat.fdr[,"Delta"]<-delta
	mat.fdr[,"p0"]<-p0
	vec.order<-as.numeric(na.exclude(vec.false[order(d)]))
	for(i in 1:le.delta){
		mat.fdr[i,"j1"]<-j1<-ifelse(any(d.diff[j0:m]>=delta[i]),
			j0-1+min(which(d.diff[j0:m]>=delta[i])),m+1)
		mat.fdr[i,"cutup"]<-ifelse(j1!=m+1,d.sort[j1],Inf)
		mat.fdr[i,"j2"]<-j2<-ifelse(any(d.diff[1:(j0-1)]<= -delta[i]) & j0!=1,
			max(which(d.diff[1:(j0-1)]<=-delta[i])),0)
		mat.fdr[i,"cutlow"]<-ifelse(j2!=0,d.sort[j2],-Inf)
		mat.fdr[i,"Called"]<-m-j1+1+j2
		mat.fdr[i,"False"]<-ifelse(j1==m+1,0,vec.order[j1])+ifelse(j2==0,0,vec.order[j2])
		mat.fdr[i,"FDR"]<-min(p0*mat.fdr[i,"False"]/max(mat.fdr[i,"Called"],1),1)
	}
	mat.fdr
}	
		 
 

