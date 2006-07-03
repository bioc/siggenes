find.a0.old<-function(data,x,y,paired=FALSE,mat.samp=NULL,B=100,balanced=FALSE,na.rm=FALSE,
		delta=0.9,alpha=(0:9)/10,include.0=TRUE,p0=NA,stable=TRUE,number.int=139,
		rand=NA,plot.legend=TRUE){
    	rs.out<-rs.cal(data,x,y,paired=paired,mat.samp=mat.samp,bal=balanced,B=B,na.rm=na.rm,
		rand=rand)
    	r<-rs.out$r        
    	s<-rs.out$s        
    	r.perm<-rs.out$r.perm
    	s.perm<-rs.out$s.perm
    	mat.samp<-rs.out$mat.samp
    	vec.a0<-quantile(s,alpha,na.rm=TRUE)
    	if(include.0)
        	vec.a0<-c(0,vec.a0)
    	n.genes<-length(na.exclude(r))     
    	sig.a0<-NULL
    	mat.post<-NULL
    	Z.norm<-qnorm(((1:n.genes)-3/8)/(n.genes+0.25))  
    	for(i in 1:length(vec.a0)){     
        	Z<-sort(r/(s+vec.a0[i]))    
        	z<-sort(as.vector(r.perm/(s.perm+vec.a0[i])))  
        	z.norm<-approx(Z,Z.norm,z,rule=2)$y  
        	mat.ratio<-ratio.est(Z.norm,z.norm,p0=p0,number.int=number.int)$mat.post 
        	mat.post<-rbind(mat.post,cbind(a0=rep(vec.a0[i],nrow(mat.ratio)),mat.ratio))  
        	sig.a0[i]<-sum(mat.ratio[which(mat.ratio[,"posterior"]>=delta),"success"]) 
    	}
    	logit.post<-log(mat.post[,"posterior"]/(1-mat.post[,"posterior"]))   
    	plot(mat.post[which(mat.post[,"a0"]==vec.a0[1]),"center"],
		logit.post[which(mat.post[,"a0"]==vec.a0[1])],
		main="Transformed Z values vs. Logit of the Posterior",xlab="transformed Z values",
        	ylab="logit(posterior)",type="l",xlim=c(-4,4),
		ylim=c(0,max(logit.post[which(logit.post!=Inf)])+0.5))
    	if(any(logit.post==Inf))
		cat("Warning: Some of the logit posterior probabilities are Inf.",
			"These probabilities are not plotted.","\n","\n")
    	for(i in 2:length(vec.a0))       
        	lines(mat.post[which(mat.post[,"a0"]==vec.a0[i]),"center"],
			logit.post[which(mat.post[,"a0"]==vec.a0[i])],col=i)
    	abline(h=log(delta/(1-delta)),lty=4)    
    	vec.a0.names<-NULL
    	vec.a0.name<-NULL
    	for(i in 1:length(alpha)){    
        	vec.a0.name[i]<-paste(c("a0=",round(vec.a0[i+1],4)," (alpha=",alpha[i],")"),
			collapse="")
        	vec.a0.names[i]<-paste(c("alpha=",alpha[i]," (",sig.a0[ifelse(include.0,i+1,i)],")"),
			collapse="")
    	}
    	if(include.0){
        	vec.a0.name<-c("a0=0",vec.a0.name)
        	vec.a0.names<-c(paste(c("a0=0 (",sig.a0[1],")"),collapse=""),vec.a0.names)
    	}
    	if(plot.legend) 
        	legend(-1.1,max(logit.post)+.6,legend=vec.a0.names,lty=1,cex=0.8,col=1:length(vec.a0),
			bty="n")
    	names(sig.a0)<-vec.a0.name
    	cat("\n","Number of significant genes for some a0:","\n")   
    	print(sig.a0)
    	a0<-vec.a0[which(sig.a0==max(sig.a0))][1]
    	cat("\n","Suggested choice for a0:",round(a0,4))
    	if(a0!=0)
        	cat("   (the",names(a0),"quantile of the s-values)")  
    	cat("\n")
    	structure(list(r=r,s=s,r.perm=r.perm,s.perm=s.perm,mat.samp=mat.samp,sig.a0=sig.a0,a0=a0,
		delta=delta,vec.a0=vec.a0,x=x,y=y))
} 
 

