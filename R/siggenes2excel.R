`siggenes2excel` <-
function(object,delta,file,excel.version=1,n.digits=3,what="both",entrez=FALSE,
		chip="",quote=FALSE){
	if(!excel.version%in%c(1,2))
		stop("'excel.version must be either 1 and 2.")
	sep<-ifelse(excel.version==1,",",";")
	dec<-ifelse(excel.version==1,".",",")
	suffix<-tolower(substring(file,nchar(file)-3,nchar(file)))
	if(suffix!=".csv"){
		file<-paste(file,"csv",sep=".")
		warning("Since the suffix of 'file' is not 'csv' this suffix is added ",
			"to 'file'.",call.=FALSE)
	}
	summary(object,delta,n.digits=n.digits,what=what,entrez=entrez,chip=chip,
		file=file,sep=sep,quote=quote,dec=dec)
}

