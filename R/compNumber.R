`compNumber` <-
function(z,post,p0,B,delta=0.9,vec.pos=NULL,vec.neg=NULL){
	if(any(delta<=0 | delta>1))
		stop("The delta values must be between 0 and 1.")
	z.sort<-sort(z)
	z.order<-order(z)
	post<-post[z.order]
	if(length(vec.pos)==0)
		probs<-1/(B*(1-post)/p0+1)
	else{
		vec.pos<-vec.pos[z.order]
		vec.neg<-vec.neg[z.order]
	}
	n.delta<-length(delta)
	m<-length(z)
	mat.delta<-matrix(0,n.delta,5)
	rownames(mat.delta)<-1:n.delta
	colnames(mat.delta)<-c("Delta","Number","FDR","CL","CU")
	mat.delta[,1]<-delta
	for(i in 1:n.delta){
		if(any(z.sort<0) & post[1]>=delta[i]){
			tmp<-post[z.sort<0]
			j1<-min(which(tmp<delta[i]))-1
			f1<-if(length(vec.neg)==0) sum(1/probs[1:j1]-1)/B
				else vec.neg[j1]
			mat.delta[i,4]<-z.sort[j1]
		}
		else{
			j1<-0
			f1<-0
			mat.delta[i,4]<- -Inf
		}
		if(post[m]>=delta[i]){
			tmp<-rev(post[z.sort>=0])
			j2<-min(which(tmp<delta[i]))-1
			f2<-if(length(vec.pos)==0) sum(1/probs[(m-j2+1):m]-1)/B
				else vec.pos[m-j2+1]
			mat.delta[i,5]<-z.sort[m-j2+1]
		}
		else{
			j2<-0
			f2<-0
			mat.delta[i,5]<-Inf
		}
		mat.delta[i,2]<-j1+j2
		mat.delta[i,3]<-min(1,p0*(f1+f2)/max(1,j1+j2))
	}
	mat.delta		
}

