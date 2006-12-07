q.value.old<-function(d,d.perm,p0){
	n.with.NA<-length(d) 
	d<-na.exclude(d)      
	d.perm<-na.exclude(d.perm)
	m<-length(d)        
	B<-length(d.perm)/length(d)
	mat.addup<-matrix(c(abs(d),abs(d.perm),rep(1,m),rep(0,m*B),sign(d),sign(d.perm)),m*(B+1),3)
	mat.addup<-mat.addup[order(mat.addup[,1]),1:3]  
	vec.false<-m-(which(mat.addup[,2]==1)-(1:m))/B  
	q.value<-p0*vec.false[1]/m  
	for(i in 2:m)
		q.value[i]<-min(p0*vec.false[i]/(m-i+1),q.value[i-1])
	mat.qvalue<-cbind(d=d[order(abs(d))],false=vec.false,q.value=q.value) 
	if(n.with.NA!=m)  
		mat.qvalue<-rbind(mat.qvalue,matrix(NA,n.with.NA-m,3))
	structure(list(mat.addup=mat.addup,vec.false=vec.false,q.value=q.value,mat.qvalue=mat.qvalue))
}
	 
 

