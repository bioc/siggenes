chisqClassSplitted<-function(data,cl,n.cat,n.split,check=TRUE,withNA=FALSE){
	n.row<-nrow(data)
	ints<-round(seq(1,n.row+1,le=n.split+1))
	ints<-unique(ints)
	stats<-numeric(n.row)
	for(i in 1:n.split)
		stats[ints[i]:(ints[i+1]-1)] <- chisqClass(data[ints[i]:(ints[i+1]-1),],
			cl,n.cat,check=check, withNA=withNA)
	stats
}