pairt.cl.transform<-function(mat,n.row){
	if(!all(dim(mat)==c(n.row,2)))
		stop("y must be either a vector or a ncol(data) x 2 matrix.")
	if(mode(as.matrix(mat))=="character"){
		vec1<-if(!all(substring(paste(mat[,1],collapse=""),1:sum(nchar(mat[,1])),
			1:sum(nchar(mat[,1])))%in%c(as.character(0:9),".","-"))) 
			as.character(mat[,1]) else as.numeric(mat[,1])
		vec2<-if(!all(substring(paste(mat[,2],collapse=""),1:sum(nchar(mat[,2])),
			1:sum(nchar(mat[,2])))%in%c(as.character(0:9),".","-"))) 
			as.character(mat[,2]) else as.numeric(mat[,2])
		mat<-data.frame(vec1,vec2)
	}
	tab1<-table(mat[,1])
	tab2<-table(mat[,2])
	n.tab1<-length(tab1)
	n.tab2<-length(tab2)
	if(n.tab1!=2 && n.tab2!=2)
		stop("One of the columns of y must contain 2 different values.")
	if(n.tab1!=n.row/2 && n.tab2!=n.row/2)
		stop("One of the columns of y must contain ",n.row/2," different values.")
	if(any(table(mat[,1],mat[,2])!=rep(1,n.row)))
		stop("There is something wrong with the y matrix.")
	vec.sign<-if(n.tab1==2) mat[,1] else mat[,2]
	vec.pair<-if(n.tab1==n.row/2) mat[,1] else mat[,2]
	sort.sign<-sort(unique(vec.sign))
	sort.pair<-sort(unique(vec.pair))
	if(!all(sort.sign==c(-1,1))){
		warning("Expected values of one of the columns of y are -1 and 1.","\n",
			as.character(sort.sign[1])," is thus set to -1, and ", 
			as.character(sort.sign[2])," to 1.",call.=FALSE)
		vec.sign<-ifelse(vec.sign==sort.sign[1],-1,1)
	}
	if(!all(sort.pair==1:(n.row/2))){
		mnc<-matrix(0,n.row/2,2)
		vec.tmp<-numeric(length(vec.pair))
                for(i in 1:(n.row/2)){
                	vec.tmp[vec.pair==sort.pair[i]]<-i
                	mnc[i,]<-c(as.character(sort.pair[i]),i)
                }
		vec.pair<-vec.tmp
		warning("Expected values of one of the columns are the integers between 1 and ",
			n.row/2,".","\n","The new class labels are thus ",paste(mnc[,2],"(was ",
			mnc[,1],ifelse(mnc[,2]!=mnc[nrow(mnc),2],"), ",")."),sep="",collapse=""),
			call.=FALSE)
	}
	y<-as.numeric(vec.sign)*as.numeric(vec.pair)
	y
}
	 
 

