na.replace.cont<-function(X,na.method="mean"){
	if(!na.method%in%c("mean","median"))
		stop("na.method must be either mean or median.")
	FUN<-match.fun(na.method)
	if(is.vector(X))
		X<-replace(X,which(is.na(X)),FUN(X,na.rm=TRUE))
	else{
		for (i in 1:nrow(X)) 
			X[i,]<-replace(X[i,],which(is.na(X[i,])),FUN(X[i,],na.rm=TRUE))
	}
	return(X)
} 
 

