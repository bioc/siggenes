args.sam<-function(method){
	if(!method%in%c("print","summary","plot","identify"))
		stop("method must be either print, summary, plot or identify.")
	if(method=="print"){
		cat("function(x,delta=NULL,n.digits=3)","\n\n")
		cat("x          a SAM object","\n",
                    "delta      a numeric value or vector specifying one or a set of Deltas.","\n",
		    "           If NULL, the Delta values of x will be used","\n",
		    "n.digits   integer specifying the number of decimal places in the output",
			sep="","\n")
	}
	if(method=="summary"){
		cat("function(object,delta=NULL,n.digits=5,file=\"\",what=\"both\")","\n\n")
		cat("object     a SAM object","\n",
		    "delta      a numeric value or vector specifying one or a set of Deltas.\n",
		    "           If NULL or a vector, summary will do the same as print with\n",
		    "           additional information. If delta is a value not only statistics\n",
		    "           as the number of significant genes and the FDR will be shown\n",
		    "           but also gene-specific information about the genes called\n",
		    "           significant when using this value of delta -- see what\n",
		    "n.digits   integer specifying the number of decimal places in the output\n",
                    "what       either \"both\", \"stats\" or \"genes\". Only available when\n",
		    "           delta is a numeric value. If \"stats\", general information is shown.\n",
                    "           If \"genes\", gene-specific information is shown.\n",
                    "           If \"both\", both general and specific information is shown\n",
		    "ll         if TRUE, both the locus links and the symbols of the genes\n",
		    "           will be added to the output\n",
                    "chip       the chip type used in this analysis. Only used, if ll=TRUE.\n",
		    "           If the argument 'data' of sam(data,cl,...) has been specified by\n",
		    "           an exprSet object, chip need not to be specified\n", 
		    "file       name of a file in which the information should be stored.\n",
		    "           By default the information is shown in the R window. Only\n", 
		    "           available when delta is a numeric value\n",
		    "sep        the field separator string used when output is stored in 'file'\n",
		    "quote      see ?write.table",
			sep="","\n","\n")
		cat("Output:","\n","If delta is NULL or a vector: number of significant genes, FDR etc.",
			"\nIf delta is a numerical value, the output will consist of \n",
			"  row.sig.genes   the rows of the data matrix containing the significant genes\n",
			"  mat.fdr         a vector containing the general information\n",
			"  mat.sig         a matrix containing gene-specific information\n",sep="")
	}
	if(method=="plot"){
		cat("function(x,y,...)","\n\n")
		cat("x    a SAM object\n",
		    "y    the delta value(s), i.e. either a numeric value or vector. If y is a\n",
                    "     numeric value, a SAM Plot for this delta value will be generated. If NULL\n",
		    "     or a vector, the Delta plots, i.e. a visualization of the table created by\n",
		    "     summary or print, are generated\n",
                    "...  further options. If y is a numeric value, see ?sam.plot2. Otherwise, see,\n",
		    "     ?delta.plot\n",sep="")
	}
	if(method=="identify"){
		cat("function(x,showText=TRUE,getInfo=TRUE,pos=4,cex=0.8,add.xy=numeric(2),\n",
			"n.digits=4,ask=TRUE,ll=TRUE,browse=FALSE,chip=\"\",...)","\n\n",sep="")
		cat("x          a SAM object\n",
		    "showText   should the gene name be plotted near the identified point?\n",
                    "getInfo    if TRUE, information of the gene corresponding to the identified\n",
		    "           point will be displayed\n",
	            "pos        the place where the gene name is plotted. Default is 4 (to the\n",
		    "           right of the point). For details see ?text\n",
                    "cex        the relative size of the plotted text. For details see ?par\n",
                    "add.xy     a vector of length 2. The text is plotted add.xy[1] units in\n",
		    "           x-direction and add.xy[2] units in the y-direction from the\n",
                    "           position where it is plotted by default\n",
                    "n.digits   integer specifying the number of decimal places in the output\n",
                    "ask        if TRUE, the user will be asked before the next point can be\n",
		    "           identified\n",
                    "ll         if TRUE, both the locus links and the symbols of the genes\n",
		    "           will be added to the output\n",        
                    "browse     if TRUE, the NCBI webpage corresponding to the locus link of\n",
		    "           the identified point is opened. Only used, if ll=TRUE\n",
                    "chip       the chip type used in this analysis. Only used, if ll=TRUE.\n",
		    "           If the argument 'data' of sam(data,cl,...) has been specified by\n",
		    "           an exprSet object, chip need not to be specified\n",
                    "...        further options such as col for the plotted text. See ?text\n",
			sep="")
	}
}   
 

