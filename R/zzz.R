.First.lib <- function(libname, pkgname) {
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
	message("\n","Latest changes in siggenes:","\n",
		"Version 1.9.1: Added a new, much faster version of sam.snp.","\n",
		"Version 1.9.9: Added a new, completely revised version of the EBAM functions.",
		"\n","\n",
		"Note: The EBAM update is not completed yet. This in particular means that","\n",
		"      the vignette has not been updated yet. Please see the help files of","\n",
		"       find.a0 and ebam for examples how to use these functions.","\n","\n",
		"Still to come: An updated vignette, EBAM analysis with Wilcoxon rank sums,","\n",
		"               EBAM analysis for categorical data","\n")
}


