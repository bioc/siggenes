check.chipname<-function(chip,obj.chip,cdfname=NULL){
		if(chip=="" & obj.chip=="" & is.null(cdfname))
			stop("The chip type must be specified either by the SAM object",
				"\n or by 'chipname' or 'cdfname'.")
		if(chip=="")
			chip<-obj.chip
		if(chip!=obj.chip & obj.chip!="")
			stop("'chipname' differs from the chip type of the SAM object.")
		if(!is.null(cdfname)){
			if(!is.character(cdfname))
				stop("'cdfname' must be a character string.")
			tmp<-tolower(cdfname)
			tmp<-gsub("_","",tmp)
			tmp<-gsub("-","",tmp)
			tmp<-gsub(" ","",tmp)
			if(chip!="" & tmp!=chip)
				stop("'chipname' differs from the chip type specified by",
					" 'cdfname'.")
			chip<-tmp
		}
		chip
}



