q.value.wilc<-function(W,p0,n.x,n.y,paired=FALSE){
    	m<-length(na.exclude(W))    
    	n<-ifelse(paired,n.x,n.x+n.y)  
    	W.mean<-ifelse(paired,n*(n+1)/4,n.x*(n+1)/2)  
    	W.min<-ifelse(paired,0,n.x*(n.x+1)/2)
    	W.max<-ifelse(paired,n*(n+1)/2,n.x*(2*n.y+n.x+1)/2)  
    	p.value<-numeric(W.max-W.min+1)  
    	nsig<-numeric(W.max-W.min+1)     
    	for(i in W.min:floor(W.mean)){  
        	nsig[i-W.min+1]<-sum(W<=i | W>=W.max+W.min-i,na.rm=TRUE)
        	nsig[length(nsig)+W.min-i]<-sum(W<=i | W>=W.max+W.min-i,na.rm=TRUE)
        	p.value[i-W.min+1]<-min(1,2*ifelse(paired,psignrank(i,n),pwilcox(i-W.min,n.x,n.y)))
        	p.value[length(p.value)+W.min-i]<-min(1,2*ifelse(paired,psignrank(i,n),
			pwilcox(i-W.min,n.x,n.y)))
    	}
    	numeric(W.max-W.min+1)->vec.q.value     
   	vec.q.value[floor(W.mean-W.min+1)]<-p0*p.value[floor(W.mean-W.min+1)]  
    	for(i in (floor(W.mean)-1):W.min)   
        	vec.q.value[i-W.min+1]<-min(p0*p.value[i-W.min+1]/(nsig[i-W.min+1]/m),
			vec.q.value[i-W.min+2])
    	vec.q.value[(ceiling(W.mean):W.max)-W.min+1]<-rev(vec.q.value[(W.min:floor(W.mean))-W.min+1])
    	q.value<-numeric(length(W))  
    	for(i in W.min:W.max)
        	q.value[which(W==i)]<-vec.q.value[i-W.min+1]
    	mat.qvalue<-cbind(W,q.value)[order(abs(W-W.mean),na.last=TRUE),]   
    	structure(list(mat.qvalue=mat.qvalue,p.value=p.value,nsig=nsig,q.value=q.value,
		vec.q.value=vec.q.value,W.min=W.min,W.max=W.max,W.mean=W.mean))
}
 
 

