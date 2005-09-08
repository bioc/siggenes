sam.sampler<-function(n.x,n.y,B,paired=FALSE,balanced=FALSE,rand=NA,file.out=NA){
    if(!is.na(rand))
        set.seed(rand)
    if(!paired){
        if(!balanced){
            mat.samp<-matrix(0,B,n.x+n.y)
            for(i in 1:B)
                mat.samp[i,]<-sample(c(rep(1,n.x),rep(0,n.y)),n.x+n.y)
        }
        if(balanced){
            mat.samp.x<-matrix(0,B,n.x)
            mat.samp.y<-matrix(0,B,n.y)
            for(i in 1:B){
                mat.samp.x[i,]<-sample(rep(c(1,0),ceiling(n.x/2)),n.x)
                mat.samp.y[i,]<-sample(rep(c(1,0),ceiling(n.y/2)),n.y)
            }
            mat.samp<-cbind(mat.samp.x,mat.samp.y)
        }
    }
    if(paired){
        if(n.x!=n.y)
            stop("x must be equal to y")
        mat.samp<-matrix(0,B,n.x)
        for(i in 1:B)
            mat.samp[i,]<-sample(c(-1,1),n.x,replace=TRUE)
    }
    if(!is.na(file.out))
        write.table(mat.samp,file=file.out,sep="\t")
    return(mat.samp)
}
    
 
 

