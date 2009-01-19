quantiles<-function(x,prob){
    if(any(prob>1 | prob<0))
        stop("Probabilities must be between 0 and 1")
    x.sort<-sort(x)
    prob.exclude<-which(prob>0 & prob<1)
    nprob<-length(x)*prob[prob.exclude]
    quan<-numeric(length(prob))
    if(any(prob==0))
        quan[which(prob==0)]<-min(x)
    if(any(prob==1))
        quan[which(prob==1)]<-max(x)
    int<-which(nprob==round(nprob))
    not.int<-which(nprob!=round(nprob))
    if(length(not.int)>0)
        quan[prob.exclude][not.int]<-x.sort[ceiling(nprob[not.int])]
    if(length(int)>0)
        quan[prob.exclude][int]<-.5*(x.sort[nprob[int]]+x.sort[nprob[int]+1])
    return(quan)
}

