samControl <- function(delta=NULL, n.delta=10, p0=NA, lambda=seq(0,0.95,0.05), ncs.value="max",
		ncs.weights=NULL, q.version=1){
	list(delta=delta, n.delta=n.delta, p0=p0, lambda=lambda, ncs.value=ncs.value,
		ncs.weights=ncs.weights, q.version=q.version)
}

ebamControl <- function(p0=NA, p0.estimation=c("splines","interval","adhoc"), lambda=NULL,
		ncs.value="max", use.weights=FALSE){
	list(p0=p0, p0.estimation=match.arg(p0.estimation), lambda=lambda, ncs.value=ncs.value,
		use.weights=use.weights)
}

find.a0Control <- function(p0.estimation=c("splines", "adhoc", "interval"), lambda=NULL,
		ncs.value="max", use.weights=FALSE, n.chunk=5, n.interval=139, 
		df.ratio=NULL){
	list(p0.estimation=match.arg(p0.estimation), lambda=lambda, ncs.value=ncs.value, 
		use.weights=use.weights, n.chunk=n.chunk, n.interval=n.interval,
		df.ratio=df.ratio)
}

 