adjust.for.mt<-function(data,cl,var.equal=FALSE,eb=FALSE,wilc=FALSE){
	if(is.data.frame(cl) || is.matrix(cl))
		cl<-pairt.cl.transform(cl,ncol(data))
	ana.type<-ifelse(eb,"EBAM","SAM")
	n.cl<-length(cl)
	if(n.cl!=ncol(data))
		stop("The length of cl must be equal to the number of columns of data.")
	lev<-unique(cl)
	uni.cl<-length(lev)
	if(any(lev<0)){
		if(any(lev==0) | length(lev[lev>0])!=length(lev[lev<0]))
			stop(paste("Negative values are only allowed for the paired analysis. Or,", "\n",
				"if a paired analysis should be done: There is something wrong with the class labels."))
		uni.cl.abs<-uni.cl/2
	}
	else
		uni.cl.abs<-uni.cl
	type.mt<-NULL
	if(uni.cl==1){
		type.mt<-"pairt"
		X<-matrix(0,nrow(data),2*n.cl)
		X[,seq(1,2*n.cl-1,2)]<-as.matrix(data)
		msg<-paste(ana.type,"Analysis for the One-Class Case",
			if(wilc) "Using Wilcoxon Signed Rank Statistics","\n\n")
		if(lev!=1)
			warning("Expected class label is 1, ",lev," is thus set to 1.",call.=FALSE)
		cl.mt<-rep(c(1,0),n.cl)
	}
	if(uni.cl==2){
		msg<-paste(ana.type,"Analysis for the Two-Class Unpaired Case",
			ifelse(!wilc,paste("Assuming",ifelse(var.equal,"Equal","Unequal"),"Variances"),
			"Using Wilcoxon Rank Sums"),"\n","\n")
		if(any(table(cl)<2))
			stop("Each group must consist of at least two samples.")
		type.mt<-ifelse(var.equal,"t.equalvar","t")
		X<-as.matrix(data)
		if(!all(sort(lev)==0:1)){
			warning("Expected class labels are 0 and 1. ",sort(lev)[1]," is set to 0, and ",
				sort(lev)[2]," is set to 1.",call.=FALSE)
			cl[which(cl==sort(lev)[1])]<-0
			cl[which(cl==sort(lev)[2])]<-1
		}
		cl.mt<-as.numeric(cl)
		
	}
	if(uni.cl==2*uni.cl.abs & uni.cl==n.cl){
		msg<-paste(ana.type,"Analysis for the Two-Class Paired Case",
			if(wilc) "Using Wilcoxon Signed Rank Scores","\n\n")
		type.mt<-"pairt"
		sort.cl<-sort(cl,index=TRUE)
		if(!all(sort.cl$x==c(-uni.cl.abs:-1,1:uni.cl.abs)))
			stop(paste("The class labels must be the integers between ",-uni.cl.abs,
				" and 1, and between 1 and ",uni.cl.abs,".",sep=""))
		x<-sort.cl$ix[(uni.cl.abs+1):uni.cl]
		y<-sort.cl$ix[uni.cl.abs:1]
		index<-as.vector(rbind(x,y))
		X<-as.matrix(data[,index])	
		cl.mt<-rep(c(1,0),n.cl/2)
	}
	if(uni.cl>2 & all(lev>=0)){
		if(wilc)
			stop("Currently no rank-based SAM version for the multi-class case available.")
		msg<-paste(ana.type,"Analysis for the Multi-Class Case with",uni.cl,"Classes","\n","\n")	
		type.mt<-"f"
		if(any(table(cl)<2))
			stop("Each group must consists of at least two samples.","\n","\n")
		sort.cl<-sort(lev)
		cl.new<-cl
		if(!all(sort.cl==1:uni.cl)){
			mnc<-matrix(0,uni.cl,2)
			for(i in 1:uni.cl){
				cl.new[which(cl==sort.cl[i])]<-i
				mnc[i,]<-c(sort.cl[i],i)
			}
			warning("Expected class labels are the integers between 1 and ",uni.cl,".",
				"\n","The new class labels are thus ",paste(mnc[,2],"(was ",mnc[,1],
				ifelse(mnc[,2]!=mnc[nrow(mnc),2],"), ",")."),sep="",collapse=""),
				call.=FALSE)
		}
		X<-as.matrix(data)
		cl.mt<-as.numeric(cl.new)-1
	}
	if(is.null(type.mt))
		stop("There is something wrong with the class labels.")
	#mode(X)<-"numeric"
	structure(list(X=X,cl.mt=cl.mt,type.mt=type.mt,msg=msg))
}
 
 

