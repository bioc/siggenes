# Copyright (c) 2002 Holger Schwender

# This program produces the SAM Plot for one or several delta. This function is constructed as helpfunction.
# So it doesn't check for incorrect use. That's why this function should carefully be used, if one likes to
# use it.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that
# 1. SAM was introduced in Tsuher, V., Tibshirani, R., and Chu, G. (2001), Significance analysis
#    of microarrays applied to the ionizing radiation response, PNAS,98, 5116-5121,
# 2. there is a patent pending for the SAM technology at Stanford University,
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# d.sort: vector of sorted observed d-values
# d.bar: vector of (sorted) expected d-values
# delta: value or vector (it is possible to choose more than one delta) for which the SAM Plot should
#        be generated.
# pty.square: if TRUE, a square SAM Plot is generated with x- and y-axes having the same range
# main: the title of the plot, which appears above the plot
# make.legend: if FALSE, a legend for deltas will be added to the SAM Plot. make.legend should only be set to TRUE,
#              if several deltas are used and neither vec.fdr nor sig.genes is specified
# vec.fdr: if vec.fdr is specified, a SAM Plot for the corresponding delta will be generated. This option
#          should only be used for one delta. 
# sig.genes: specifying sig.genes will only make sense, if vec.fdr is specified. If sig.genes is specified,
#            the points in the SAM Plot which correspond to the values of sig.genes will be marked with green
#            color
# wilc: if TRUE, a SAM Plot for sam.wilc() will be generated. If FALSE, a SAM Plot for sam() is made.
# use.numbers: if TRUE, the symbols for the points will be the numbers of observations which correspond to this
#              points. Otherwise circles are used. This option should only be used, if wilc=TRUE. 
# count: a vector with the numbers of observations which correspond to the points in the plot. This should
#        only be used, if wilc=TRUE.


sam.plotter<-function(d.sort,d.bar,delta,pty.square=TRUE,main="SAM Plot",color=1,make.legend=FALSE,vec.fdr=NULL,
            sig.genes=NULL,wilc=FALSE,use.numbers=FALSE,count=NULL){
    lim.min.x<-min(d.bar-.1,na.rm=TRUE)  # some limits are set for the plot
    lim.max.x<-max(d.bar+.1,na.rm=TRUE)
    lim.min.y<-min(d.sort-.1,na.rm=TRUE)
    lim.max.y<-max(d.sort+.1,na.rm=TRUE)
    if(pty.square){
        par(pty="s",lab=c(10,10,7))
        lim.min.x<-min(lim.min.x,lim.min.y)  # other limits will be needed, if a
        lim.min.y<-min(lim.min.x,lim.min.y)  # square SAM Plot is desired
        lim.max.x<-max(lim.max.x,lim.max.y)
        lim.max.y<-max(lim.max.x,lim.max.y)
    }
    lab.names<-ifelse(wilc,"W values","d(i)")  # make some dinstinction between the labels
    symbol<-ifelse(use.numbers,"n","p")        # for sam() and sam.wilc()
    plot(d.bar,d.sort,main=main,xlab=paste(c("expected",lab.names),collapse=" "),ylab=paste(c("observed",lab.names),collapse=" "),
        type=symbol,xlim=c(lim.min.x,lim.max.x),ylim=c(lim.min.y,lim.max.y))  # points in the SAM Plot
    if(use.numbers)    # numbers can be used as symbols in sam.wilc()
        text(d.bar,d.sort,count,cex=.8)
    abline(0,1)
    for(i in 1:length(delta)){   # the delta lines
        abline(delta[i],1,lty=2,col=color[i])
        abline(-delta[i],1,lty=2,col=color[i])
    }
    if(make.legend)   # a legend is made
        legend(lim.min.x,lim.max.y,legend=c("delta",as.character(delta)),lty=2,col=c(0,color),bty="n",cex=.8)
    if(!is.null(vec.fdr)){        # the cutlines are plotted
        abline(h=vec.fdr[6],lty=5)  # cutlow
        abline(h=vec.fdr[7],lty=5)  # cutlup
        # text is added to the plot
        text(rep(lim.min.x,6),seq(lim.max.y,lim.max.y-(lim.max.y-lim.min.y)/4,length=6),c("cutlow:","cutup:",
            "p0:","significant:","false:","FDR:"),adj=0,cex=.75)
        text(rep(lim.min.x+(lim.max.x-lim.min.x)/ifelse(pty.square,4.8,6),6),seq(lim.max.y,lim.max.y-(lim.max.y-lim.min.y)/4,length=6),
            round(vec.fdr[c(6,7,2,4,3,5)],3),adj=0,cex=.75)
        if(use.numbers)
            text(d.bar[sig.genes],d.sort[sig.genes],count[sig.genes],cex=.8,col=3)
        else
            points(d.bar[sig.genes],d.sort[sig.genes],col=3)
    }
}
