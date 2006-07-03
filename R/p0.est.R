p0.est<-function(d,d.perm,lambda=1,vec.lambda=(0:95)/100){
	if(lambda>1 || lambda<0)   
		stop("lambda has to be between 0 and 1.")
	if(lambda!=1)       
		vec.lambda<-lambda
	vec.p0<-NULL
	d<-na.exclude(d)
	m<-length(d)
	quan<-quantiles(na.exclude(d.perm),c(vec.lambda/2,1-rev(vec.lambda)/2))  
	for(i in 1:length(vec.lambda))
		vec.p0[i]<-sum(d>quan[i] & d<quan[length(quan)-i+1])/((1-vec.lambda[i])*m)  
	if(lambda!=1){
		p0<-min(vec.p0,1)      
		return(list(p0=p0,vec.p0=vec.p0))
	}
	spline.out<-smooth.spline(vec.lambda,vec.p0,w=1-vec.lambda,df=3)   
	p0<-min(predict(spline.out,1)$y,1)   
	structure(list(p0=p0,spline.out=spline.out,vec.p0=vec.p0))
}
