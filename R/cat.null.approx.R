"cat.null.approx" <-
function(d,n.cl,n.cat){
	df<-(n.cl-1)*(n.cat-1)
	n.var<-length(d)
	d.bar<-qchisq(((1:n.var)-0.5)/n.var,df)
	p.value<-pchisq(d,df,lower.tail=FALSE)
	vec.false<-p.value*n.var
	return(list(d.bar=d.bar,p.value=p.value,vec.false=vec.false))
}



	

