build.dperm<-function(X,tmp.samp,type.mt,s0,n.row,le.cl){
	B<-nrow(tmp.samp)
	max.cl<-max(tmp.samp[1,])+1
	d.perm<-matrix(0,n.row,B)
	for(i in 1:B){
		tmp<-.C("get_stat_num_denum", as.double(X), as.integer(n.row), 
        	as.integer(le.cl), as.integer(tmp.samp[i,]), as.double(.mt.naNUM), 
        	t.num = double(n.row), t.denum=double(n.row),as.character(c(type.mt,"abs","y")), 
		as.integer(max.cl), PACKAGE = "multtest")
		d.perm[,i]<-sort(tmp$t.num/(tmp$t.denum+s0),na.last=TRUE)
	}
	d.perm
}	 
 

