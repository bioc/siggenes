"cat.null" <-
function(data,mat.samp,d,n.subset,n.cat){
	B<-nrow(mat.samp)
	vecB<-c(seq(1,B,n.subset),B+1)
	vecB<-unique(vecB)
	n.B<-length(vecB)-1
	n.var<-nrow(data)
	vec.dperm<-vec.false<-numeric(n.var)
	d.rank<-rank(-d,ties.method="first")
	for(i in 1:n.B){
		tmp<-compPermStat(data,mat.samp[vecB[i]:(vecB[i+1]-1),,drop=FALSE],n.cat)
		tmp<-apply(tmp,2,sort)
		vec.dperm<-vec.dperm+rowSums(tmp)
		tmp2<-c(as.vector(tmp),d)
		vec.false<-vec.false+rank(-tmp2,ties.method="first")[length(tmp)+(1:n.var)]-d.rank
	}
	d.bar<-vec.dperm/B
	vec.false<-vec.false/B
	p.value<-vec.false/n.var
	return(list(d.bar=d.bar,vec.false=vec.false,p.value=p.value))
}

