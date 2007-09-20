.onLoad <- function(libname, pkgname){
	require(methods)
}

.onAttach <- function(libname, pkgname) {
	message("\n","Latest major changes in siggenes:","\n",
		"Version 1.9.27: EBAM update completed. Added new vignette.","\n",
		"Version 1.11.7: Improved version of EBAM for categorical data.","\n",
		"Version 1.11.9: Added findDelta for determining the number of genes called\n",
                "                differentially expressed for a given FDR, and vice versa.", 
		"\n")
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
}


