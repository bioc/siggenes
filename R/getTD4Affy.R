getTD4Affy<-function(genenames,cdfname,isMap=FALSE){
	tmp<-paste("https://www.affymetrix.com/analysis/netaffx/",if(isMap) "mapping",
		"fullrecord.affx?pk=",cdfname,":",genenames,sep="")
	td<-paste("<TD><A HREF=\"",tmp,"\">",genenames,"</A></TD>",sep="")
	td
}



