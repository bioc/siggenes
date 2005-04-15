make.plotname<-function(device="png"){
	now<-Sys.time()
	plotname<-paste("samplot",paste(substring(now,c(3,6,9),c(4,7,10)),collapse=""),
		paste(sample(0:9,6),collapse=""),sep="_")
	plotname<-paste(plotname,device,sep=".")
	plotname
}



