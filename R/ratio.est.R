ratio.est<-function(Z.norm,z.norm,p0=NA,stable=TRUE,number.int=139){
    	require(splines)  
    	min.int<-floor(100*min(Z.norm))/100     
    	max.int<-ceiling(100*max(Z.norm))/100   
    	interval<-seq(min.int,max.int,length=number.int+1)    
    	center<-(interval[2]-interval[1])/2+interval[-length(interval)]
    	bin.Z<-cut(Z.norm,interval,include.lowest=TRUE)       
   	bin.z<-cut(z.norm,interval,include.lowest=TRUE)       
    	success<-tabulate(bin.Z,length(levels(bin.Z)))    
    	failure<-tabulate(bin.z,length(levels(bin.z)))    
    	n<-success+failure    
    	p<-success/n          
    	log.bino<-lgamma(n+1)-(lgamma(success+1)+lgamma(n-success+1)) 
    	ns.out<-ns(center,5)        
    	mat.repeat<-na.exclude(cbind(log.bino,ns.out,n,success,p,center))  
    	mat.repeat<-as.data.frame(mat.repeat)                              
	names(mat.repeat)<-c("log.bino","x1","x2","x3","x4","x5","n","success","p","center")
      	attach(mat.repeat) 
    	optim.out<-optim(rep(0,6),neglogLik.repeat,method="BFGS")
    	b<-as.vector(optim.out$par)      
    	mat.model<-as.matrix(cbind(1,mat.repeat[,2:6])) 
    	pi.Z<-exp(mat.model%*%b)/(1+exp(mat.model%*%b))    
    	B<-length(z.norm)/length(Z.norm)     
    	if(is.na(p0)){
        	if(stable)
            		p0<-p0.est(Z.norm,z.norm)$p0
        	else
            		p0<-min((B*pi.Z)/(1-pi.Z))
        }
    	posterior<-1-p0*(1-pi.Z)/(B*pi.Z) 
    	posterior[which(posterior<0)]<-0 
    	mat.post<-cbind(mat.repeat[,c("center","success")],posterior)
    	structure(list(mat.repeat=mat.repeat,optim.out=optim.out,p0=p0,mat.post=mat.post))
}
 
 

