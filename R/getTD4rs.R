getTD4rs<-function(snps,refsnp){
	if(is.null(names(refsnp)))
		stop("The vector (or data frame) refsnp has no (row)names.")
	if(any(!snps%in%names(refsnp)))
		stop("Some of the SNPs specified by snps are missing in refsnp.")
	refs<-refsnp[snps]
	ids.nolink<-which(is.na(refs) | substring(refs,1,2)!="rs")
	tmp<-paste("http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=",refs,sep="")
	td<-paste("<TD><A HREF=\"",tmp,"\">",refs,"</A></TD>",sep="")
	if(length(ids.nolink)>0)
		td[ids.nolink]<-"<TD>&nbsp;</TD>"
	td
}

