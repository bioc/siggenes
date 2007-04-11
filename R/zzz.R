.OnLoad <- function(libname, pkgname){
	require(methods)
}

.onAttach <- function(libname, pkgname) {
	message("\n","Latest changes in siggenes:","\n",
		"Version 1.9.1: Added a new, much faster version of sam.snp.","\n",
		"Version 1.9.9: Added a new, completely revised version of the EBAM functions.","\n",
		"Version 1.9.27: EBAM update completed. Added new vignette.",
		"\n")
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
}


