"compPermStat" <-
function(X,P,n.cat){
	if(missing(n.cat))
		n.cat<-max(X)
	n.class<-max(P)
	n.obs<-ncol(X)
	vecX<-vector("list",n.cat)
	vecP<-vector("list",n.class)
	for(i in 1:n.cat)
		vecX[[i]]<-X==i
	for(i in 1:n.class)
		vecP[[i]]<-P==i
	vecXs<-lapply(vecX,rowSums)
	vecRs<-table(P[1,])
	#vecS<-vector("list",n.cat*n.class)
	matStat<-matrix(0,nrow(X),nrow(P))
	if(n.class>2){
		for(i in 1:n.class){
			for(j in 1:n.cat){
				tmp<-vecX[[j]]%*%t(vecP[[i]])
				#tmp<-tmp*tmp
				tmp2<-vecXs[[j]]*vecRs[i]/n.obs
				#vecS[[(i-1)*n.cat+j]]<-tmp*tmp/tmp2
				matStat<-matStat+tmp*tmp/tmp2	
			}
		}
	}
	else{
		for(j in 1:n.cat){
			tmp<-vecX[[j]]%*%t(vecP[[1]])
			tmp2<-vecXs[[j]]*vecRs[1]/n.obs
			matStat<-matStat+tmp*tmp/tmp2
			tmp<-vecXs[[j]]-tmp
			tmp2<-vecXs[[j]]-tmp2
			matStat<-matStat+tmp*tmp/tmp2
		}
	}
	#matStat<-matrix(0,nrow(X),nrow(P))
	#for(i in 1:length(vecS))
	#	matStat<-matStat+vecS[[i]]
	matStat-n.obs
}

