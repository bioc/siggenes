.onLoad <- function(libname, pkgname){
	require(methods)
}

.onAttach <- function(libname, pkgname) {
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui"){
        	addVigs2WinMenu("siggenes")
		winMenuAddItem("Vignettes/siggenes","Rnews Paper",
			"openPDF(file.path(.find.package('siggenes'),'doc','siggenesRnews.pdf'))")
	}
}


