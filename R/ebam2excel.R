`ebam2excel` <-
function(object,delta,file,excel.version=1,n.digits=4,what="both",entrez=FALSE,
		chip="",quote=FALSE){
	if(!is(object,"EBAM"))
		stop("object must be an object of class EBAM.")
	siggenes2excel(object,delta,file,excel.version=excel.version,n.digits=n.digits,
		what=what,entrez=entrez,chip=chip,quote=quote)
}

