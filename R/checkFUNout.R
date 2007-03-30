`checkFUNout` <-
function(fun.out){
	tmp.names<-names(fun.out)
	if(any(!c("z","ratio")%in%tmp.names))
		stop("The output of the function specified by method must contain objects\n",
			"called z and ratio.")
	if(length(fun.out$z)!=length(fun.out$ratio))
		stop("z and ratio must have the same length.")
	exact<-ifelse(any(tmp.names=="vec.pos"),TRUE,FALSE)
	a0<-if(any(tmp.names=="a0")) fun.out$a0 else numeric(0)
	msg<-if(any(tmp.names=="msg")) fun.out$msg else ""
	if(any(tmp.names=="mat.samp")){
		mat.samp<-fun.out$mat.samp
		B<-nrow(mat.samp)
	}
	else{
		mat.samp<-matrix(numeric(0))
		B<-NA
	}
	if(!exact)
		vec.pos<-vec.neg<-numeric(0)
	else{
		vec.pos<-fun.out$vec.pos
		if(length(vec.pos)!=length(fun.out$z))
			stop("vec.pos must have the same length as z.")
		if(any(tmp.names=="vec.neg")){
			vec.neg<-fun.out$vec.neg
			if(length(vec.neg)!=length(vec.pos))
				stop("vec.neg must have the same length as vec.pos.")
		}
		else{
			twosided<-ifelse(any(fun.out$z<0),TRUE,FALSE)
			vec.neg<-if(twosided) vec.pos else numeric(length(vec.pos))
			warning("vec.neg is not specified even though vec.pos is specified.\n",
				"Since ",ifelse(twosided,"some","none")," of the z values ",
				ifelse(twosided,"are","is")," smaller than zero, ",
				"vec.neg is set to ",ifelse(twosided,"vec.pos","0"),".",
				call.=FALSE)
		}
	}
	return(list(vec.pos=vec.pos,vec.neg=vec.neg,exact=exact,B=B,mat.samp=mat.samp,a0=a0,msg=msg))
}

