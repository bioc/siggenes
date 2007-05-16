wilc.ebam<-function(data,cl,approx50=TRUE,ties.method=c("min","random","max"),use.offset=TRUE,
		df.glm=5,use.row=FALSE,rand=NA){
	require(splines)
	if(!is.na(rand))
		set.seed(rand)
	data<-as.matrix(data)
	mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("No missing values allowed.")
	adjust.out<-adjust.for.mt(data,cl,wilc=TRUE,eb=TRUE)
	X<-adjust.out$X
	cl.mt<-adjust.out$cl.mt
	type.mt<-adjust.out$type.mt
	msg<-adjust.out$msg
	rm(data,adjust.out)
	n.row<-nrow(X)
	n.cl<-length(cl.mt)
	ties.method<-match.arg(ties.method)
	if(type.mt=="t"){
		n0<-sum(cl.mt==0)
		n1<-sum(cl.mt==1)
		W.mean<-n1*(n.cl+1)/2
		W.min<-n1*(n1+1)/2
		W.max<-n1*(2*n0+n1+1)/2
		W.var<-n1*n0*(n.cl+1)/12
		if(use.row)
			W<-rowWilcoxon(X,cl.mt)
		else{
			X.rank<-matrix(0,n.row,n.cl)
			for(i in 1:n.row)
				X.rank[i,]<-rank(X[i,],ties.method=ties.method)
			W<-rowSums(X.rank[,cl.mt==1])
		}
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n0<50 & n1<50)
			approx50<-FALSE
		if(!approx50){
			W.null<-dwilcox(W-W.min,n1,n0)
			p.value<-2*pwilcox(W.mean-abs(W-W.mean)-W.min,n1,n0)
			p.value[p.value>1]<-1
			tmp<-sort(unique(W-W.min))
			f.null<-dwilcox(tmp,n1,n0)
		}

		else{
			W.null<-dnorm(W.norm)
			p.value<-2*(1-pnorm(abs(W.norm)))
			f.null<-dnorm(sort(unique(W.norm)))
		}
	}
	if(type.mt=="pairt"){
		n.cl<-n.cl/2
		X<-X[,2*(1:n.cl)-1]-X[,2*(1:n.cl)]
		W.max<-n.cl*(n.cl+1)/2
		W.mean<-W.max/2
		W.var<-n.cl*(n.cl+1)*(2*n.cl+1)/24
		if(sum(X==0)>0){
			warning("There are ",sum(X==0)," observation pairs having a difference of zero.",
				"\n","These differences are randomly set to either 1e-06 or -1e-06.",
				call.=FALSE)
			X[X==0]<-sample(c(1e-06,-1e-06),sum(X==0),rep=TRUE)
		}
		if(use.row)
			W<-rowWilcoxon(X,rep(1,n.cl))
		else{
			X.rank<-matrix(0,n.row,n.cl)
			for(i in 1:n.row)
				X.rank[i,]<-rank(abs(X[i,]),ties.method=ties.method)
			W<-rowSums(X.rank*(X>0))
		}
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n.cl<50)
			approx50<-FALSE
		if(!approx50){
			W.null<-dsignrank(W,n.cl)
			p.value<-2+psignrank(W.mean-abs(W-W.mean),n.cl)	
			p.value[p.value>1]<-1
			f.null<-dsignrank(sort(unique(W)),n.cl)
		}
		else{
			W.null<-dnorm(W.norm)
			p.value<-2*(1-pnorm(abs(W.norm)))
			f.null<-dnorm(sort(unique(W.norm)))
		}
	}
	vec.pos<-n.row*p.value/2
	tabW<-table(W)
	valW<-as.numeric(names(tabW))
	offset.value<-if(use.offset) log(f.null) else rep(0,length(valW))
	glm.out<-glm(tabW~ns(valW,df.glm)+offset(offset.value),family=poisson)
	f.W<-glm.out$fitted/n.row
	W.fitted<-f.W[as.character(W)]
	list(z=W.norm,ratio=W.null/W.fitted,vec.pos=vec.pos,vec.neg=vec.pos,msg=msg)
}
	
	
		




	