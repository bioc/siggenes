getQuantiles<-function(n.knots,mode){
	den<-(n.knots+1)/2
	tmp1<-(1:floor(n.knots/2))*mode/den
	tmp2<-1-(1-mode)*(ceiling(n.knots/2):1)/den
	c(tmp1,tmp2)
}

