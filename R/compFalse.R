`compFalse` <-
function(z,z.perm,z.rank,n.row){
	ids<-which(z.perm<0)
	n.neg<-length(ids)
	if(n.neg>0)
		vec.neg<-rank(c(z.perm[ids],-abs(z)),ties.method="first")[n.neg+(1:n.row)]-z.rank
	else
		vec.neg<-numeric(n.row)
	ids<-which(z.perm>=0)
	vec.pos<-rank(-c(z.perm[ids],abs(z)),ties.method="first")[length(ids)+(1:n.row)]-z.rank
	return(list(vec.pos=vec.pos,vec.neg=vec.neg))
}

