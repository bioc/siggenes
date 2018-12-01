nclass.wand<-function(x,level=1){
	# requireNamespace("KernSmooth", quietly=TRUE)
	ceiling(diff(range(x))/KernSmooth::dpih(x,level=level))
}


