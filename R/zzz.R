.First.lib <- function(libname, pkgname) {
	require(Biobase) || stop("cannot load siggenes without Biobase")
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
	cat("\n","Changes in siggenes:","\n",
		"Version 1.2.7: Added sam2excel and sam2html","\n",
		"Version 1.2.17: Added help.sam (help files for SAM-methods plot,",
		" summary, ...)","\n",
		"Version 1.5.3: Now exprSet objects can also be used in the EBAM functions","\n",
		"               for specifying 'data' and 'cl'.","\n",
		"For all changes previous to 1.2.7, see the siggenes vignette.",
		"\n",sep="")

}
