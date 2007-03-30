`truncZ` <-
function(z.perm,z.min,z.max){
	z.perm[z.perm<z.min]<-z.min
	z.perm[z.perm>z.max]<-z.max
	z.perm
}

