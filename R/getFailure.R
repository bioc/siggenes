`getFailure` <-
function(z.perm,z,interval,z.norm=NULL,z.range=NULL,n.interval=139){
	if(is.null(z.norm))
		z.perm<-truncZ(z.perm,z.range[1],z.range[2])
	else{
		z.sort<-sort(z)
		z.perm<-approx(z.sort,z.norm,z.perm,rule=2)$y
	}
	bin<-cut(z.perm,interval,include.lowest=TRUE)
	tabulate(bin,n.interval)
}

