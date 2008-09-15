buildSNPannotation<-function(pkg,rs=TRUE,allele=TRUE,gene=TRUE,chromosome=FALSE,
		position=FALSE,strand=FALSE,cytoband=FALSE,lib.loc=NULL,others=NULL,
		subset=NULL,pattern=NULL,na.rm=TRUE){
	require(pkg,character.only=TRUE,lib.loc=lib.loc) || stop(paste("Package",pkg,
		"not available."))
	conn<-db(get(pkg))
	what<-c("man_fsetid","dbsnp_rs_id","chrom","physical_pos","strand","cytoband",
		"allele_a","allele_b","gene_assoc")
	cn<-c("Probe-Set-ID","RefSNP","Chromosome","Position","Strand","Cytoband",
		"Allele_A","Allele_B","Gene")
	interest<-c(TRUE,rs,chromosome,position,strand,cytoband,allele,allele,gene)
	what<-what[interest]
	cn<-cn[interest]
	if(!is.null(others)){
		what<-unique(c(what,others))
		cn<-unique(c(cn,others))
	}
	if(any(!what%in%dbListFields(conn,"featureSet")))
		stop("Some of the specified annotations seem to be not available.")
	sql<-paste("SELECT", paste(what,collapse=", "), "FROM featureSet")
	if(!is.null(pattern))
		sql<-paste(sql," WHERE man_fsetid LIKE '",pattern,"'",sep="")
	out<-dbGetQuery(conn,sql)
	rn<-out[,1]
	rownames(out)<-rn
	colnames(out)<-cn
	str<-out$Strand
	if(any(str%in%(0:1))){
		str[str==0]<-"+"
		str[str==1]<-"-"
		out$Strand<-str
	}
	if(!is.null(subset)){
		ids<-match(subset,rn)
		if(any(is.na(ids))){
			warning(sum(is.na(ids))," of the ",length(ids)," Probe-Set-IDs specified",
				" by 'subset'\n","are not available in the annotation package.")
			if(na.rm){
				subset<-subset[!is.na(ids)]
		 		ids<-ids[!is.na(ids)]
			}
		}
		out<-out[ids,]
		rownames(out)<-subset
	}
	out
}
	