`sam2excel` <-
function(object,delta,file,excel.version=1,n.digits=3,what="both",entrez=FALSE,
		chip="",quote=FALSE){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	siggenes2excel(object,delta,file,excel.version=excel.version,n.digits=n.digits,
		what=what,entrez=entrez,chip=chip,quote=quote)
}

