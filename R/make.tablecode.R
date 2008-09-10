make.tablecode<-function(genenames,entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,
		chipname="",cdfname=NULL,dataframe=NULL,new.window=TRUE,tableborder=1,
		name1stcol="Name",refsnp=NULL,which.refseq="NM",load=TRUE){
	require(annotate)
	if(!is.null(refsnp))
		entrez<-refseq<-symbol<-omim<-ug<-FALSE
	if(any(c(entrez,refseq,symbol,omim,ug))){
		if(chipname=="" & is.null(cdfname))
			stop("Either 'chipname' or 'cdfname' must be specified.")
	}
	if(!is.null(cdfname) && is.null(refsnp)){
		require(affy)
		clean<-cleancdfname(cdfname,addcdf=FALSE)
		if(chipname=="")
			chipname<-clean
		if(clean!=chipname)
			stop("'chipname' and 'cdfname' do not specify the same chip.")
		tmp<-new("AffyBatch",cdfName=cdfname,annotation=clean)
		gN<-geneNames(tmp)
		if(!all(genenames%in%gN))
			stop("Some of the 'genenames' do not specify genes of the ",
				cdfname," chip.")
		tr<-getTD4Affy(genenames,cdfname)
	}
	else
		tr<-paste("<TD>",genenames,"</TD>",sep="")
	tr<-paste(tr,"\n")
	th<-name1stcol
	if(symbol){
		sym<-unlist(lookUp(genenames,chipname,"SYMBOL",load=load))
		sym[is.na(sym)]<-"&nbsp;"
		sym.cols<-paste("<TD>",sym,"</TD>\n",sep="")
		th<-c(th,"Symbol")
		tr<-paste(tr,sym.cols)
	}
	if(!is.null(refsnp)){
		rs<-getTD4rs(genenames,refsnp)
		th<-c(th,"RefSNP")
		tr<-paste(tr,rs,"\n")
	}
	if(entrez){
		LL<-unlist(lookUp(genenames,chipname,"ENTREZID",load=load))
		LL[is.na(LL)]<-"&nbsp;"
		LL.cols<-annotate:::getTDRows(LL,"en")
		th<-c(th,"LocusLink")
		tr<-paste(tr,LL.cols,"\n")
	}
	if(refseq){
		tmp.refseq<-lookUp(genenames,chipname,"REFSEQ",load=load)
		which.refseq<-c(tolower(which.refseq),toupper(which.refseq))
		nm.select<-function(x) x[substring(x,1,2)%in%which.refseq]
		RefSeq<-lapply(tmp.refseq,nm.select)
		RefSeq[lapply(RefSeq,length)==0]<-"&nbsp;"
		refseq.cols<-annotate:::getTDRows(RefSeq,"GB")
		th<-c(th,"RefSeq")
		tr<-paste(tr,refseq.cols,"\n")
	}
	any.na<-function(x){any(is.na(x))}
	if(omim){
		OMIM<-lookUp(genenames,chipname,"OMIM",load=load)
		OMIM[unlist(lapply(OMIM,any.na))]<-"&nbsp;"
		OMIM.cols<-annotate:::getTDRows(OMIM,"omim")
		th<-c(th,"OMIM")
		tr<-paste(tr,OMIM.cols,"\n")
	}
	if(ug){
		UG<-lookUp(genenames,chipname,"UNIGENE",load=load)
		UG[lapply(UG,length)==0]<-"&nbsp;"
		# Notloesung:
		UG<-unlist(lapply(UG,function(x) x[1]))
		UG.cols<-annotate:::getTDRows(UG,"ug")
		th<-c(th,"UniGene")
		tr<-paste(tr,UG.cols,"\n")
	}
	if(new.window)
		tr<-add.target2href(tr)
	if(!is.null(dataframe)){
		if(!is.data.frame(dataframe))
			stop("'dataframe' must be a data.frame.")
		if(nrow(dataframe)!=length(genenames))
			stop("'The number of rows of 'data.frame' differ from the length",
				" of 'genenames'.")
		th<-c(th,colnames(dataframe))
		for(i in 1:ncol(dataframe)){
			tmp<-paste("<TD align=\"right\">",dataframe[,i],"</TD>",sep="")
			tr<-paste(tr,tmp,"\n")
		}
	}
	tr<-paste("<TR>\n",tr,"</TR>\n")
	trth<-paste(c("<TR>\n",paste("<TH>",th,"</TH>\n",sep=""),"</TR>"),collapse=" ")
	tr<-c(trth,"\n",tr)
	tr<-c(paste("<table border=",tableborder," align=center cellpadding=5px>\n",
		sep=""),tr,"</table>\n")
	tr	
}




