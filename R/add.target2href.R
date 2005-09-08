add.target2href<-function(tr,split="\">"){
	tmp<-strsplit(tr,split)
	for(i in 1:length(tmp)){
		tmp1<-tmp[[i]]
		if(length(tmp1)>1)
			tr[i]<-paste(tmp1[1],paste("\" target=\"_blank\">",tmp1[-1],
				sep="",collapse=""),sep="")
	}
	tr
}



