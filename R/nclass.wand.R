nclass.wand<-function(x,level=1){
	require(KernSmooth)
	ceiling(diff(range(x))/dpih(x,level=level))
}


