fudge2<-function(r,s,alpha=seq(0,1,0.05),include.zero=TRUE){
	if (max(alpha) > 1 || min(alpha) < 0) 
        	stop("alpha has to be between 0 and 1")
    	if (any(round(100 * alpha, 10) != round(100 * alpha, 0))) {
        	warning("At least one alpha is not a percentile. Only the first two decimal digits",
			 " are retained.")
        alpha <- signif(alpha, 2)
    	}
	if(length(alpha)==1){
		s.zero<-quantile(s,alpha)
		msg<-paste("s0 =",round(s.zero,4)," (The",100*alpha,"% quantile of the s values.) \n \n")
		invisible(return(list(s.zero=s.zero,alpha.hat=alpha,vec.cv=NULL,msg=msg)))
	}
	fudge.quan<-quantile(s,alpha)
	if(include.zero)
		fudge.quan<-c(0,fudge.quan)
	n.alpha<-length(fudge.quan)
	d.mat<-r/outer(s,fudge.quan,"+")
	n.uni.s <- length(unique(s))
	if(n.uni.s<25)
		stop("For the computation of the fugde factor,","\n",
			"there should be at least 25 genes with differing standard deviations.")
    	n.int <- ifelse(n.uni.s > 500, 101, floor(n.uni.s/5))
    	quan <- quantile(s, seq(0, 1, le = n.int))
	quan<-unique(round(quan,8))
	n.int<-length(quan)
	int.s<-as.numeric(cut(s,quan,include.lowest=TRUE,right=FALSE))
	mad.mat<-matrix(0,n.int-1,ncol(d.mat))
	for(i in 1:(n.int-1)){
		mad.mat[i,]<-apply(d.mat[which(int.s==i),,drop=FALSE],2,mad)
		#med.s<-which(d.mat[which(int.s==i),1]==median(d.mat[which(int.s==i),1]))
		#med.mat[which(int.s==i),]<-d.mat[med.s,]
	}
	cv<-function(x){
		sd(x)/mean(x)
	}
	vec.cv<-apply(mad.mat,2,cv)
	which.min<-which(vec.cv==min(vec.cv))
	if(include.zero & which.min==1){
		msg<-"s0 = 0 \n \n"
		s.zero <- 0
		invisible(return(list(s.zero=s.zero,vec.cv=vec.cv,msg=msg)))
	}
	s.zero<-fudge.quan[which.min]
	if(include.zero)
		which.min<-which.min-1
	alpha.hat<-alpha[which.min]
        msg<-paste("s0 =", round(s.zero, 4), " (The", 100 * alpha.hat, 
            "% quantile of the s values.)", "\n", "\n")
	invisible(return(list(alpha.hat=alpha.hat,s.zero=s.zero,vec.cv=vec.cv,msg=msg)))
}



 
 

