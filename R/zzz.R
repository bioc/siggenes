.OnLoad <- function(libname, pkgname){
	require(methods)
}

.onAttach <- function(libname, pkgname) {
	message("\n","Latest changes in siggenes:","\n",
		"Version 1.9.1: Added a new, much faster version of sam.snp.","\n",
		"Version 1.9.9: Added a new, completely revised version of the EBAM functions.","\n",
		"Version 1.9.22: Added NAMESPACE to siggenes.",
		"\n","\n",
		"Note: Although the EBAM update is completed, the vignette has not been updated","\n",
		"      yet. Please see the help files of find.a0 and ebam for examples of","\n",
		"      how to use these functions.","\n")
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
}


