neglogLik.repeat<-function(b){
	sum(-success*(b[1]+b[2]*x1+b[3]*x2+b[4]*x3+b[5]*x4+b[6]*x5)
		+n*log(1+exp(b[1]+b[2]*x1+b[3]*x2+b[4]*x3+b[5]*x4+b[6]*x5)))
} 
 

