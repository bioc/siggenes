link.genes<-function(genenames,filename,entrez=TRUE,refseq=TRUE,symbol=TRUE,omim=FALSE,
		ug=FALSE,fullname=FALSE,which.refseq="NM",chipname="",cdfname=NULL,
		dataframe=NULL,title=NULL,bg.col="white",text.col="black",link.col="blue",
		tableborder=1,new.window=TRUE,load=TRUE){
	tr<-make.tablecode(genenames,entrez=entrez,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		fullname=fullname,chipname=chipname,cdfname=cdfname,which.refseq=which.refseq,
		dataframe=dataframe,tableborder=tableborder,new.window=new.window,load=load)
	suffix<-tolower(substring(filename,nchar(filename)-4,nchar(filename)))
	if(suffix!=".html"){
		filename<-paste(filename,"html",sep=".")
		warning("Since the suffix of 'filename' is not 'html' '.html' is added",
			" to 'filename'.",call.=FALSE)
	}
	bg.col<-col2hex(bg.col)
	text.col<-col2hex(text.col)
	link.col<-col2hex(link.col)
	if(is.null(title))
		title<-"Links for a Set of Genes to Public Repositories"
	outfile<-file(filename,"w")
	cat("<html>","<head>","<title>Links to Public Repositories</title>","</head>",
		paste("<body bgcolor=",bg.col," text=",text.col," link=",link.col,">",
		sep=""),
		paste("<h1 align=center>",title,"</h1>",sep=""),
		paste("<style type=text/css>",sep=""),
		"p{ margin-top: 2px; margin-bottom: 2px;}",
		"</style>",
		tr,"</body>","</html>",sep="\n",file=outfile)
	close(outfile)
}




