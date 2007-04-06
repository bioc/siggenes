`getSuccesses` <-
function(z,n.interval=139){
	min.int<-floor(100*min(z))/100
	max.int<-ceiling(100*max(z))/100
	interval<-seq(min.int,max.int,length=n.interval+1)
	center<-(interval[2]-interval[1])/2 + interval[-length(interval)]
	bin<-cut(z,interval,include.lowest=TRUE)
	success<-tabulate(bin,n.interval)
	return(list(success=success,interval=interval,center=center))
}

