R.fold.old<-function(data,x,y,na.rm=FALSE){
    X<-as.matrix(data[,c(x,y)])
    mode(X)<-"numeric"
    mean.x <- rowMeans(X[,1:length(x)],na.rm=na.rm)  
    mean.y <- rowMeans(X[,(length(x)+1):ncol(X)], na.rm=na.rm)  
    vec.R.fold<-mean.x/mean.y  # compute the fold changes
    vec.R.fold[which(mean.x<=0 | mean.y<=0)]<-NA   
    mat.R.fold<-cbind(mean.x=mean.x,mean.y=mean.y,R.fold=vec.R.fold) 
    mat.R.fold
}
 
 

