filterALL<-function(){
	fp<-.find.package("siggenes")
	fp<-file.path(fp,"doc/filterALL2.R")
	source(fp)
	filterALL2()
}


