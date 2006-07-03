na.replace.dist<-function(x){
	x[is.na(x)]<-sample(na.exclude(sort(unique(x))),sum(is.na(x)),replace=TRUE,prob=table(x))
	x
} 
 

