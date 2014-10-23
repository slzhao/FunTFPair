##' convertId
##' 
##' A function to convert ID based on the biomaRt package.
##' 
##' A function to convert ID based on the biomaRt package..
##' 
##' @param keepMultipleId Logical. Indicate keep the multiple target IDs related to one source ID or not.
##' @param keepNoId Logical. Indicate keep the source IDs without target IDs or not.
##' @param verbose Logical. Indicate report extra information on progress or not.
##' @inheritParams biomaRt::getBM
##' @inheritParams biomaRt::useMart
##' @inheritParams newIdMatrix
##' @importFrom biomaRt getBM useMart
##' @export
##' @examples temp<-cbind(rnorm(10),rnorm(10))
##' row.names(temp)<-c("Q04837","P0C0L4","P0C0L5","O75379","Q13068","A2MYD1","P60709","P30462","P30475","P30479")
##' colnames(temp)<-c("Exp1","Exp2")
##' convertId(temp,filters="uniprot_swissprot_accession",keepMultipleId=TRUE)
##' \dontrun{
##' temp<-cbind(rnorm(5000),rnorm(5000),rnorm(5000),rnorm(5000),rnorm(5000),rnorm(5000))
##' row.names(temp)<-1000:5999
##' colnames(temp)<-c("Control1","Control2","Control3","Treatment1","Treatment2","Treatment3")
##' convertId(temp,filters="entrezgene",attributes =c("entrezgene","uniprot_swissprot_accession"),keepNoId=FALSE)
##' }
convertId<-function(x,dataset="hsapiens_gene_ensembl",filters="uniprot_swissprot_accession",attributes =c(filters,"entrezgene"),genesKept=c('foldchange','first','random','var','abs'),keepNoId=T,keepMultipleId=F,verbose=F) {
#	if (! require("biomaRt")) {
#		cat("biomaRt package is needed but not installed in this computer. Will install it from bioconductor.\n")
#		flush.console()
#		source("http://bioconductor.org/biocLite.R")
#		biocLite("biomaRt")
#		if (!require(biomaRt)) {stop("Package biomaRt can't be installed")}
#	}
	if (missing(genesKept)) {
		genesKept<-"var"
	} else {
		genesKept<-match.arg(genesKept)
	}
	if (verbose) {
		cat("Now conectting with ensembl. Internet acess is needed and it may use 30 seconds.\n")
		flush.console()
	}
	ensembl = useMart("ensembl",dataset=dataset)
	
	oldId<-row.names(x)
	newIdTable<-getBM(attributes =attributes,filters=filters,values=oldId,mart = ensembl)
	newIdTable<-newIdTable[which(newIdTable[,1]!="" & newIdTable[,2]!=""),]
	
	temp1<-which(oldId %in% newIdTable[,1])
	temp2<-nrow(x)-length(temp1)
	xNoId<-NULL
	if (keepNoId) {
		if (verbose) {
			cat(paste("No ID, Keep: ",temp2," genes can't find their ",attributes[2]," ID. They will be attched at the end of data with their original ID.\n",sep=""))
		}
		xNoId<-x[-temp1,,drop=F]
	} else {
		if (verbose) {
			cat(paste("No ID, Discard: ",temp2," genes can't find their ",attributes[2]," ID. They will be discard.\n",sep=""))
		}
	}
	
	temp2<-split(newIdTable,newIdTable[,1])
	if (keepMultipleId) {
		newIdTable<-sapply(temp2,function(x) return(paste(x[,2],collapse=";")))
		if (verbose) {
			cat(paste("Multiple IDs, Keep: ",length(which(sapply(temp2,function(x) nrow(x))>=2))," genes have more than one ",attributes[2]," IDs. All of these ",attributes[2]," IDs will be stored.\n",sep=""))
		}
	} else {
		newIdTable<-sapply(temp2,function(x) return(x[1,2]))
		if (verbose) {
			cat(paste("Multiple IDs, Discard: ",length(which(sapply(temp2,function(x) nrow(x))>=2))," genes have more than one ",attributes[2]," IDs. Only the first ",attributes[2]," ID for each gene will be used.\n",sep=""))
		}
	}
	names(newIdTable)<-names(temp2)
	
	result<-newIdMatrix(x,genesKept=genesKept,convertIdTable=newIdTable)
	
	if (keepMultipleId) {
		temp<-strsplit(row.names(result),";")
		temp1<-sapply(temp,length)
		temp2<-unlist(temp)
		result<-result[rep(1:length(temp1),temp1),]
		row.names(result)<-temp2
	}
	if (keepNoId) {
		result<-rbind(result,xNoId)
	}
	return(result)
}

##' newIdMatrix
##' 
##' A function to convert ID.
##' 
##' A function to convert ID.
##' 
##' @param x the expression data matrix.
##' @param convertIdTable A vector. The names should be the source IDs, and the values should be the target IDs.
##' @param genesKept The method to select target gene in more than one targets. "var"/"foldchange"/"abs" means selecting the gene with largest variation/fold change/absolute value. "first" means selecting the first target and "random" means randomly selection.
##' @export
##' @examples convertIdTable<-paste("New",c(1,2,2,2,1,3,4,4,5,5))
##' names(convertIdTable)<-paste("Old",1:length(convertIdTable))
##' temp<-matrix(rnorm(20),ncol=2)
##' row.names(temp)<-names(convertIdTable)
##' colnames(temp)<-c("Exp1","Exp2")
##' newIdMatrix(temp,genesKept="foldchange",convertIdTable)
newIdMatrix<-function(x,convertIdTable,genesKept=c("var","foldchange","abs","first","random")) {
	convertIdTable<-convertIdTable[which(convertIdTable!="" & names(convertIdTable)!="")]
	x<-x[names(convertIdTable),,drop=F]
	if (missing(genesKept)) {
		genesKept<-"var"
	} else {
		genesKept<-match.arg(genesKept)
	}
	if (genesKept=="foldchange") {
		temp<-apply(x,1,range,na.rm=T)
		testStat<-temp[2,]-temp[1,]
	} else if (genesKept=="first") {
		testStat<-rep(1,nrow(x))
		names(testStat)<-row.names(x)
	} else if (genesKept=="random") {
		testStat<-sample(1:nrow(x),nrow(x))
		names(testStat)<-row.names(x)
	} else if (genesKept=="var") {
		testStat<-apply(x,1,var,na.rm=T)
	} else if (genesKept=="abs") {
		testStat<-apply(x,1,function(y) max(abs(y),na.rm=TRUE))
	}
	temp<-split(testStat,convertIdTable)
	result<-x[unlist(sapply(temp, function(y) names(which.max(y)))),]
	row.names(result)<-names(temp)
	return(result)
}



