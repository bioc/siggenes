col2hex<-function(col){
	rgb<-col2rgb(col)[,1]
	hex<-c(0:9,LETTERS[1:6])
	paste("#",paste(hex[rgb%/%16+1],hex[rgb%%16+1],sep="",collapse=""),sep="")
}



