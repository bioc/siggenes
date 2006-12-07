.First.lib <- function(libname, pkgname) {
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
	message("\n","Changes in siggenes:","\n",
		"Version 1.5.3: Now exprSet objects can also be used in the EBAM functions","\n",
		"               for specifying 'data' and 'cl'.","\n",
                "Version 1.7.1: siggenes can also handle ExpressionSet objects.","\n",
		"Version 1.9.1: Added a new, much faster version of sam.snp.","\n")
}


