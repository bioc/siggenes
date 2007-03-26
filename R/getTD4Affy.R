getTD4Affy<-function(genenames,cdfname){
	tmp<-paste("https://www.affymetrix.com/analysis/netaffx/fullrecord.affx?pk=",
		cdfname,":",genenames,sep="")
	td<-paste("<TD><A HREF=\"",tmp,"\">",genenames,"</A></TD>",sep="")
	td
}



