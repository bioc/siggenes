# Copyright (c) 2002 Holger Schwender

# This function produces a SAM Plot for a delta, and it stores some statistics (p0, #significant genes,
# #falsely called genes, FDR) and a table of the significant genes with some statistics (d-value, s-value)
# in a file. A previous analysis with sam() must be done to use this function, but it is possible to choose
# a delta which was not used in this previous analysis.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# sam.out: the object to which the output of a previous analyis with sam() was assigned
# delta: the value of delta with which the significant genes and the FDR should be identified. This delta
#        must not necessarily be used in the previous analysis   
# data: the used data set; this data set must be the same data set as in sam(), but it could be, e.g.,
#        the unnormalized version of the data set used in sam() if this data set was normalized
# q.values: if TRUE, for each gene its q-value is computed 
# R.fold: if TRUE, the fold change of each significant gene is calculated and added to the output
# na.rm: if na.rm=FALSE, the d-values of genes with one or more missing values will be set to NA. If na.rm=TRUE, the
#        missing values will be removed during the computation of the d-values.
# pty.square: if TRUE, a square SAM Plot is generated with x- and y-axes having the same range
# file.out: output (for details see description of this function) is stored in this file. To prevent this,
#           set file.out=NA
# col.accession: if col.accession is a positive integer, this column of data is interpreted as the accession
#                number of the gene and it is added to the output. To avoid this, set col.accession=NA
# col.gene.name: if col.gene.name is a positive integer, this column of data is interpreted as the name of
#                the gene and is added to the output. To avoid this, set col.gene.name=NA
# use.numbers: if TRUE, the number of observations which correspond to a point in the SAM Plot will be used
#              as symbol for this point (instead of a circle)
# rand: if specified, the set.seed() which is used in the computation of the q-values to remove 
#       tied Wilcoxon Rank statistics is set to rand


sam.plot.old<-function(sam.out,delta,data,q.values=TRUE,R.fold=TRUE,na.rm=FALSE,pty.square=TRUE,file.out=NA,
        col.accession=NA,col.gene.name=NA,use.numbers=sam.out$use.numbers,rand=sam.out$rand){
    wilc<-ifelse(is.null(sam.out$d),TRUE,FALSE)
    if(wilc){
        d<-sam.out$W
        d.sort<-sam.out$W.sort
        d.bar<-sam.out$W.exp
        d.perm<-sam.out$table.W.exp
    }
    else{
        d<-sam.out$d
        d.bar<-sam.out$d.bar
        d.sort<-sam.out$d.sort
        d.perm<-sam.out$d.perm
    }
    p0<-sam.out$p0
    if(!any(sam.out$FDR[,1]==delta))  # checks if delta was used in previous analysis
        # if not an analysis with sam.fdr() is done, i.e. the interesting statistics are computed
        sam.fdr(d.sort,d.bar,d.perm,p0,delta=delta,med=sam.out$med.fdr,graphic=FALSE,wilc=wilc)$mat.fdr->vec.fdr
    if(any(sam.out$FDR[,1]==delta))    # if delta was used, we use the corresponding statistics saved in sam.out
        vec.fdr<-sam.out$FDR[which(sam.out$FDR[,1]==delta),]
    sig.genes<-NULL
    if(vec.fdr[9]!=0)
        sig.genes<-c(sig.genes,1:vec.fdr[9])
    if(vec.fdr[8]!=length(na.exclude(d))+1)
        sig.genes<-c(sig.genes,vec.fdr[8]:length(na.exclude(d)))
    index<-(1:length(d))                 # make an index vector of the row numbers of all genes                         
    row.sig.genes<-index[order(d)][sig.genes]   # the row numbers of the significant genes   
    if(!wilc){
        sam.plotter(d.sort,d.bar,delta,pty.square=pty.square,main=paste(c("SAM Plot for delta =",delta),collapse=""),
            vec.fdr=vec.fdr,sig.genes=sig.genes) # SAM Plot
        # the output is made
        output<-cbind("d(i)"=round(d.sort[sig.genes],4),"s(i)"=round(sam.out$s[order(d)][sig.genes],4))
        if(q.values){         #q-values are computed
            mat.qvalue<-q.value.cal(d,d.perm,p0)$mat.qvalue
            q.value<-mat.qvalue[order(mat.qvalue[,1]),3]   # sort the q.values
            output<-cbind(output,"q-value"=round(q.value[sig.genes],5))
    }}
    else{
        mat.count<-sam.out$mat.count
        fdr.for.plot<-sam.fdr(mat.count[,2],mat.count[,1],d.perm,p0,delta=delta,graphic=FALSE,wilc=TRUE)$mat.fdr
        sig.genes.plot<-c(0:fdr.for.plot[,9],fdr.for.plot[,8]:length(na.exclude(d)))
        sam.plotter(mat.count[,2],mat.count[,1],delta,pty.square=pty.square,vec.fdr=vec.fdr,sig.genes=sig.genes.plot,
            main=paste(c("SAM Plot using Wilcoxon Rank Statistics \n and delta =",delta),collapse=" "),wilc=TRUE,
            use.numbers=use.numbers,count=mat.count[,3])
        if(q.values)
            q.value<-q.value.wilc(d,p0,length(sam.out$x),length(sam.out$y),paired=sam.out$paired)$q.value
        output<-cbind("W"=d[row.sig.genes],"q-value"=round(q.value[row.sig.genes],5))
    }   
    if(!is.na(col.accession))      
        output<-cbind("access"=data[row.sig.genes,col.accession],output) 
    if(R.fold){
        fold.change<-R.fold.cal(data[row.sig.genes,],sam.out$x,sam.out$y,na.rm=na.rm)
        output<-cbind(output,"R.fold"=round(fold.change[,3],4))
    }
    if(!is.na(col.gene.name))
        output<-cbind(output,"gene"=substring(data[row.sig.genes,col.gene.name],1,50))
    output<-cbind("ID"=row.sig.genes,output)  # add the index to the output; also for an nicer output
    if(wilc)
        output<-output[nrow(output):1,]
    sam.output<-as.data.frame(output)
    if(!is.na(file.out)){  # output is stored in a file
        which.sam<-ifelse(wilc,"SAM-Wilc","SAM")
        cat("Results of",which.sam,"using delta =",round(delta,4),"\n","\n","\n","cutlow:",round(vec.fdr[6],3),"\n","cutup:",
            round(vec.fdr[7],3),"\n","p0:",round(p0,4),"\n","significant genes:",vec.fdr[4],"\n", 
            "falsely called genes:",round(vec.fdr[3],4),"\n","FDR:",round(vec.fdr[5],4),"\n","\n","\n",
            "Genes called significant",ifelse(wilc,":",paste(c("(with s0 =",round(sam.out$s0,3),"):"),collapse=" ")),
            "\n","\n",file=file.out)
	write.table(t(dimnames(sam.output)[[2]]),file=file.out,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
        write.table(sam.output,file=file.out,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
        cat("Output is stored in",file.out,"\n")
    }
    else
	print(sam.output)
    par(pty="m")    
    structure(list(vec.fdr=vec.fdr,sam.output=sam.output,row.sig.genes=row.sig.genes))
}
