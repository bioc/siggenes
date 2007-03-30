filterALL<-function(){
	pdat<-pData(ALL)
	subset<-intersect(grep("^B",as.character(pdat$BT)),
		which(pdat$mol %in% c("BCR/ABL","NEG")))
	eset<-ALL[,subset]
	require(genefilter)
	f1<-pOverA(0.25,log2(100))
	f2<-function(x) IQR(x)>0.5
	selected<-genefilter(eset,filterfun(f1,f2))
	esetSub<-eset[selected,]
	pdat<-pData(esetSub)
	esetSub$mol.biol<-as.character(esetSub$mol.biol)
	#newPData<-new("phenoData",pData=pdat,
	#	varLabels=phenoData(esetSub)@varLabels)
	#phenoData(esetSub)<-newPData
	esetSub
}