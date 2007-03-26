`computeRS` <-
function(data,cl,type){
	n.row<-nrow(data)
	le.cl<-length(cl)
	max.cl<-max(cl)+1
	tmp<-.C("get_stat_num_denum",as.double(data),as.integer(n.row),as.integer(le.cl),
		as.integer(cl),as.double(.mt.naNUM),t.num=double(n.row),t.denum=double(n.row),
		as.character(c(type,"abs","y")),as.integer(max.cl),PACKAGE="multtest")
	return(list(r=tmp$t.num,s=tmp$t.denum))
}

