#' @import Biobase
#' @name pairs2TargetRemDup
#' @title transcription factors pairs and their targets from ENCODE data.
#' @description This data set is a data.frame with 9 columns. The first 6 columns were from the sTable 3 of reference. The last 3 columns were targets for TF1 and TF2, TF1 only, TF2 only.
#' @docType data
#' @usage pairs2TargetRemDup
#' @source Gerstein MB, et al. (2012) Architecture of the human regulatory 
#' network derived from ENCODE data, Nature, Sep 6;489(7414):91-100.
NULL
##' performLimma
##' 
##' A function to perform limma analysis in prepared geo data.
##' 
##' A function to perform limma analysis in prepared geo data.
##' 
##' @param x the expression data matrix, can be result from \code{\link{prepareGeoData}}.
##' @param groups indicate the samples belong to which group.
##' @param contrasts indicate groups to be compared.
##' @param LogC logical. Indicate if a log2 transformation will be performed. If not specified, the function will try to detect it. 
##' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
##' @export
##' @examples \dontrun{
##' dataMatrix<-prepareGeoData(GEO="GDS2213")
##' limmaResult<-performLimma(dataMatrix,groups=c(rep("untreated",4),rep("aza",4),rep("TSA",4),rep("azaAndTSA",4)),contrasts=c("aza-untreated","TSA-untreated","azaAndTSA-untreated"))
##' }
performLimma<-function(x,groups,contrasts,LogC) {
#	require(limma)
	
	groups<-make.names(groups)
	x<-autoLog2(x,LogC=LogC)
	
#	eset = new("ExpressionSet", exprs=as.matrix(x))
	eset <- ExpressionSet(assayData=as.matrix(x))
	
	fl <- as.factor(groups)
	design <- model.matrix(~ fl + 0)
	colnames(design) <- levels(fl)
	fit <- lmFit(eset, design)
	cont.matrix <- makeContrasts(contrasts=contrasts, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2, 0.01)
	result<-list()
	for (i in 1:length(contrasts)) {
		result[[i]]<-topTable(fit2, adjust.method="fdr", number=99999,coef=i)
	}
	names(result)<-contrasts
	return(result)
}

##' prepareGeoData
##' 
##' A function to download GDS or GSE data and convert ID if needed.
##' 
##' A function to download GDS or GSE data and convert ID if needed.
##' 
##' @param GEO the GEO ID.
##' @param destdir indicate the path to download the files.
##' @param geneSymColumnName indicate the column name for gene symbol, which will be used in ID converter For GSE data.
##' @importFrom GEOquery getGEO Columns Table
##' @export
##' @examples dataMatrix<-prepareGeoData("GDS2213")
prepareGeoData<-function(GEO=NULL,destdir=getwd(),geneSymColumnName="Gene Symbol") {
	if (length(grep("GDS",GEO)>0)) { #gds
		gset <- getGEO(GEO=GEO,destdir=destdir)
		gsetTable<-Table(gset)
		gsetTable[,1]<-as.character(gsetTable[,1])
		gsetTable[,2]<-as.character(gsetTable[,2])
		gsetTable[,-(1:2)] <- sapply(gsetTable[,-(1:2)], as.numeric)
		#remove Control or no gene probes
		temp<-which(gsetTable[,2]=="--Control" | gsetTable[,2]=="")
		if (length(temp)) {
			gsetTable<-gsetTable[-temp,]
		}
		#remove all NA probes
		temp<-apply(gsetTable[,-(1:2)],1,function(x) length(which(is.na(x))))
		if (length(which(temp==(ncol(gsetTable)-2)))>0) {
			gsetTable<-gsetTable[-(which(temp==(ncol(gsetTable)-2))),]
		}
		temp1<-gsetTable[,-(1:2)]
		row.names(temp1)<-gsetTable[,1]
		temp2<-gsetTable[,2]
		names(temp2)<-gsetTable[,1]
		gsetTableGene<-newIdMatrix(temp1,temp2)
		pkgInfo("The experiment design:")
		print(Columns(gset))
	} else { #gse
		gset <- getGEO(GEO, GSEMatrix =TRUE,destdir =getwd())
		if (length(gset)>1) {
			pkgInfo(length(gset)," datasets were found in ",GEO,", only the first one will be used.")
		}
		gset<-gset[[1]]
		temp1<-exprs(gset)
		if (geneSymColumnName %in% colnames(fData(gset))) {
			temp2<-as.character(fData(gset)[,geneSymColumnName])
			names(temp2)<-as.character(fData(gset)[,c("ID")])
		} else {
			pkgInfo("Can't find '",geneSymColumnName,"' column for gene symbol converter. Please set geneSymColumnName parameter from the following column names:")
			print(fData(gset)[sample(1:nrow(fData(gset)),5),])
			stop()
		}

		#remove all NA probes
		temp<-apply(temp1,1,function(x) length(which(is.na(x))))
		if (length(which(temp==(ncol(temp1))))>0) {
			temp1<-temp1[-(which(temp==(ncol(temp1)))),]
			temp2<-temp2[-(which(temp==(ncol(temp1))))]
		}
		gsetTableGene<-newIdMatrix(temp1,temp2)
		pkgInfo("The experiment design:")
		print(pData(gset))
	}
	return(gsetTableGene)
}


autoLog2<-function(x,LogC) {
	if (missing(LogC)) {
		qx <- as.numeric(quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
		LogC <- (qx[5] > 100) ||
				(qx[6]-qx[1] > 50 && qx[2] > 0) ||
				(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	}
	if (LogC) {
		cat("Log2 will be performed.\n")
		x[which(x <= 0)] <- NaN
		x <- log2(x) 
	}
	return(x)
}

##' differentialAnalysis
##' 
##' A wrapper for function performLimma and enrichTarget.
##' 
##' A function to perform differential analysis for expression data, to see if the shared targets of two tf were significantly changed between two conditions.
##' 
##' @inheritParams performLimma
##' @inheritParams enrichTarget
##' @export
##' @examples \dontrun{
##' dataMatrix<-prepareGeoData(GEO="GDS2213")
##' #use TF pairs in HEPG2 cell only 
##' tfPairsUsed<-pairs2TargetRemDup[which(pairs2TargetRemDup[,3]=="HEPG2"),]
##' enrichTargetResult<-differentialAnalysis(dataMatrix,groups=c(rep("untreated",4),rep("aza",4),rep("TSA",4),rep("azaAndTSA",4)),contrasts=c("aza-untreated","TSA-untreated","azaAndTSA-untreated"),pairs2Target=tfPairsUsed)
##' enrichTargetResultNetwork1<-networkVis(enrichTargetResult[[1]])
##' enrichTargetResultNetwork3<-networkVis(enrichTargetResult[[3]])
##' }
differentialAnalysis<-function(x,groups,contrasts,LogC=NULL,pairs2Target=pairs2TargetRemDup,pCut=0.05,fCut=log2(1.5)) {
	limmaResult<-performLimma(x,groups=groups,contrasts=contrasts)
	enrichTargetResult<-enrichTarget(limmaResult,pairs2Target=pairs2Target,pCut=pCut,fCut=fCut)
	return(enrichTargetResult)
}

##' enrichTarget
##' 
##' A function to perform enrichment analysis in limma result, to see if the shared targets of two tf were enriched in significant different genes.
##' 
##' A function to perform enrichment analysis in limma result, to see if the shared targets of two tf were enriched in significant different genes.
##' 
##' @param limmaResult the result of \code{\link{performLimma}}.
##' @param pairs2Target the targets of tf pairs
##' @param pCut the cutoff of p value for significant genes
##' @param fCut the cutoff of fold change for significant genes
##' @export
##' @examples \dontrun{
##' dataMatrix<-prepareGeoData(GEO="GDS2213")
##' limmaResult<-performLimma(dataMatrix,groups=c(rep("untreated",4),rep("aza",4),rep("TSA",4),rep("azaAndTSA",4)),contrasts=c("aza-untreated","TSA-untreated","azaAndTSA-untreated"))
##' tfPairsUsed<-pairs2TargetRemDup[which(pairs2TargetRemDup[,3]=="HEPG2"),]
##' enrichTargetResult<-enrichTarget(limmaResult,pairs2Target=tfPairsUsed)
##' }
enrichTarget<-function(limmaResult,pairs2Target=pairs2TargetRemDup,pCut=0.05,fCut=log2(1.5)) {
	enrichResultList<-list()
	cat(paste0("Threre are ",length(limmaResult)," comparasions in limmaresult:\n"))
	for (i in 1:length(limmaResult)) {
		allGenes<-toupper(row.names(limmaResult[[i]]))
		sigGenes<-toupper(allGenes[which(abs(limmaResult[[i]][,1])>=fCut & limmaResult[[i]][,5]<=pCut)])
		cat(paste0("There are ",length(allGenes)," genes in comparasion ",i," and ",length(sigGenes), " of them were significant different based on log2 fold change larger than ",round(fCut,2)," and adjusted p value less than ",pCut,"\n"))
		
		resultEnrich<-matrix(NA,ncol=15,nrow=nrow(pairs2Target))
		colnames(resultEnrich)<-c("shareTargets.in.Expression","shareTargets.in.Significant","tf1OnlyTargets.in.Expression","tf1OnlyTargets.in.Significant","tf2OnlyTargets.in.Expression","tf1Targets.in.Significant","tf2OnlyTargets.in.Expression","tf1Targets.in.Significant","tf2Targets.in.Expression","tf2Targets.in.Significant","p.shareTargets","p.tf1OnlyTargets","p.tf2OnlyTargets","p.tf1Targets","p.tf2Targets")
		temp<-apply(pairs2Target[,1:3],1,function(x) paste(x,collapse="_"))
		row.names(resultEnrich)<-temp
		
		pb <- txtProgressBar(min = 1, max = nrow(pairs2Target), style = 3)
		for (x in 1:nrow(pairs2Target)) {
			setTxtProgressBar(pb, x)	
			nUsed<-NULL
			pvalue<-NULL
			for (z in c("shareTargets","tf1OnlyTargets","tf2OnlyTargets")) {
				targets<-strsplit(pairs2Target[x,z]," ")[[1]]
				mPhyper<-length(intersect(targets,allGenes))
				qPhyper<-length(intersect(targets,sigGenes))
				nPhyper<-length(allGenes)-mPhyper
				kPhyper<-length(sigGenes)
				nUsed<-c(nUsed,mPhyper,qPhyper)
				pvalue<-c(pvalue,phyper(qPhyper, mPhyper, nPhyper, kPhyper, lower.tail = F))
			}
			for (y in c(3,5)) {
				mPhyper<-nUsed[1]+nUsed[y]
				qPhyper<-nUsed[2]+nUsed[y+1]
				nPhyper<-length(allGenes)-mPhyper
				kPhyper<-length(sigGenes)
				nUsed<-c(nUsed,mPhyper,qPhyper)
				pvalue<-c(pvalue,phyper(qPhyper, mPhyper, nPhyper, kPhyper, lower.tail = F))
			}
			resultEnrich[x,]<-c(nUsed,pvalue)
		}
		cat("\n")
		resultEnrich[,11]<-p.adjust(resultEnrich[,11],method="BH")
		resultEnrich[,12]<-p.adjust(resultEnrich[,12],method="BH")
		resultEnrich[,13]<-p.adjust(resultEnrich[,13],method="BH")
		resultEnrich[,14]<-p.adjust(resultEnrich[,14],method="BH")
		resultEnrich[,15]<-p.adjust(resultEnrich[,15],method="BH")
		enrichResultList[[i]]<-resultEnrich
	}
	return(enrichResultList)
}

##' correlationAnalysis
##' 
##' A function to perform correlation analysis for expression data, to see if the shared targets of two tf have significantly higher correlations.
##' 
##' A function to perform correlation analysis for expression data, to see if the shared targets of two tf have significantly higher correlations.
##' 
##' @param expression the expression matrix, can be the result of \code{\link{prepareGeoData}}.
##' @param pairs2Target the targets of tf pairs
##' @param n the number of genes selected for prior comparison, to make the analysis faster
##' @export
##' @examples \dontrun{
##' dataMatrix3<-prepareGeoData(GEO="GSE54698",geneSymColumnName="GENE_SYMBOL")
##' dataMatrix3Used1<-dataMatrix3[,1:12]
##' dataMatrix3Used2<-dataMatrix3[,c(13:24)]
##' tfPairsUsed<-pairs2TargetRemDup[which(pairs2TargetRemDup[,3]=="HELAS3"),]
##' correlationResult3New1<-correlationAnalysis(dataMatrix3Used1,pairs2Target=tfPairsUsed,n=100)
##' correlationResult3New2<-correlationAnalysis(dataMatrix3Used2,pairs2Target=tfPairsUsed,n=100)
##' correlationResultNew1Network<-networkVis(correlationResult3New1)
##' correlationResultNew2Network<-networkVis(correlationResult3New2)
##' }
correlationAnalysis<-function(expression,pairs2Target=pairs2TargetRemDup,n=300) {
	shareTargets<-strsplit(pairs2Target[,7]," ")
	tf1OnlyTargets<-strsplit(pairs2Target[,8]," ")
	tf2OnlyTargets<-strsplit(pairs2Target[,9]," ")
	
	temp1<-(unique(unlist(shareTargets)))
	temp2<-(unique(unlist(tf1OnlyTargets)))
	temp3<-(unique(unlist(tf2OnlyTargets)))
	temp4<-intersect(row.names(expression),c(temp1,temp2,temp3))
	if (length(temp4)<=10) {
		stop(paste0("[FunTFPair] There are only ",length(temp4)," genes overlap between TF targets and expression dataset. Please check the rownames of expression dataset\n"))
	}
	expressionTarget<-expression[temp4,]
	expressionTargetCor<-cor(t(expressionTarget))
	resultCor<-matrix(NA,ncol=6,nrow=nrow(pairs2Target))
	row.names(resultCor)<-apply(pairs2Target[,1:3],1,function(x) paste(x,collapse="_"))
	colnames(resultCor)<-c("numberOfShareTargets","numberOfTf1OnlyTargets","numberOfTf2OnlyTargets","shareVsTf1OnlyPValue","shareVsTf2OnlyPValue","shareVsBackground")
		
	tempStart<-1
	tempEnd<-100
	pb <- txtProgressBar(min = tempStart, max = tempEnd, style = 3)
	
	for (x in tempStart:tempEnd) {
		setTxtProgressBar(pb, x)
		temp1<-shareTargets[x][[1]]
		temp2<-tf1OnlyTargets[x][[1]]
		temp3<-tf2OnlyTargets[x][[1]]
		sharedTargetOne<-intersect(temp1,row.names(expressionTarget))
		tf1OnlyTargetOne<-intersect(temp2,row.names(expressionTarget))
		tf2OnlyTargetOne<-intersect(temp3,row.names(expressionTarget))
		backgroundGeneOne<-sample(setdiff(row.names(expressionTarget),c(temp1,temp2,temp3)),1000)
		
		if (length(sharedTargetOne)<=5 || length(tf1OnlyTargetOne)<=5 || length(tf2OnlyTargetOne)<=5) {
			pvalue12<-NA
			pvalue13<-NA
			pvalue14<-NA
		} else {
			temp1<-targetCor(expressionTargetCor,sharedTargetOne,n)
			temp2<-targetCor(expressionTargetCor,tf1OnlyTargetOne,n)
			temp3<-targetCor(expressionTargetCor,tf2OnlyTargetOne,n)
			temp4<-targetCor(expressionTargetCor,backgroundGeneOne,n)
			pvalue12<-t.test(temp1,temp2,alternative ="g")$p.value
			pvalue13<-t.test(temp1,temp3,alternative ="g")$p.value
			pvalue14<-t.test(temp1,temp4,alternative ="g")$p.value
			if (pvalue12<=0.05 | pvalue13<=0.05 | pvalue14<=0.05) {
				temp1<-targetCor(expressionTargetCor,sharedTargetOne,99999)
				temp2<-targetCor(expressionTargetCor,tf1OnlyTargetOne,99999)
				temp3<-targetCor(expressionTargetCor,tf2OnlyTargetOne,99999)
				temp4<-targetCor(expressionTargetCor,backgroundGeneOne,99999)
				pvalue12<-t.test(temp1,temp2,alternative ="g")$p.value
				pvalue13<-t.test(temp1,temp3,alternative ="g")$p.value
				pvalue14<-t.test(temp1,temp4,alternative ="g")$p.value
			}
		}
		resultCor[x,]<-c(length(sharedTargetOne),length(tf1OnlyTargetOne),length(tf2OnlyTargetOne),pvalue12,pvalue13,pvalue14)
	}
	cat("\n")
	resultCor[,4]<-p.adjust(resultCor[,4],method="BH")
	resultCor[,5]<-p.adjust(resultCor[,5],method="BH")
	resultCor[,6]<-p.adjust(resultCor[,6],method="BH")
	return(resultCor)
}

##' networkVis
##' 
##' A function to perform functional TF network visualization
##' 
##' A function to perform functional TF network visualization
##' 
##' @param result the result of TF pairs identification, can be the result of \code{\link{differentialAnalysis}} or \code{\link{correlationAnalysis}}.
##' @param pCut the p value cutoff for genes used in network visualization
##' @import igraph
##' @export
##' @examples \dontrun{
##' dataMatrix<-prepareGeoData(GEO="GDS2213")
##' #use TF pairs in HEPG2 cell only 
##' tfPairsUsed<-pairs2TargetRemDup[which(pairs2TargetRemDup[,3]=="HEPG2"),]
##' enrichTargetResult<-differentialAnalysis(dataMatrix,groups=c(rep("untreated",4),rep("aza",4),rep("TSA",4),rep("azaAndTSA",4)),contrasts=c("aza-untreated","TSA-untreated","azaAndTSA-untreated"),pairs2Target=tfPairsUsed)
##' #Show the network for contrast 1
##' enrichTargetResultNetwork1<-networkVis(enrichTargetResult[[1]])
##' #Show the network for contrast 3, which is different from contrast 1
##' enrichTargetResultNetwork3<-networkVis(enrichTargetResult[[3]])
##' }
networkVis<-function(result,pCut=0.05) {
	if (ncol(result)>10) { #differentialAnalysis result
		resultUsed<-as.data.frame(result[,c(11,12,13)])
	} else { #correlationAnalysis result
		resultUsed<-as.data.frame(result[,c(6,4,5)])
		colnames(resultUsed)[1]<-"p.shareTargets"
	}
	#0 to min value
	resultUsed[resultUsed==0]<-min(resultUsed[resultUsed!=0],na.rm=T)
	#delete all NA pair
	resultUsed<-resultUsed[which(apply(resultUsed,1,function(x) length(which(is.na(x))))!=3),]
	
	temp<-t(sapply(strsplit(row.names(resultUsed),"_"),function(x) x))
	resultUsed<-cbind(temp,resultUsed,stringsAsFactors=F)
	
	if (ncol(result)>10) { #differentialAnalysis result
		resultUsed<-resultUsed[which(resultUsed[,4]<=pCut),]
		temp1<-sapply(split(resultUsed,as.character(resultUsed[,1]),drop = FALSE),function(x) max(x[,5]))
		temp2<-sapply(split(resultUsed,as.character(resultUsed[,2]),drop = FALSE),function(x) max(x[,6]))
		nodesP<-sapply(split(c(temp1,temp2),names(c(temp1,temp2))),max)
		nodesP<-nodesP[nodesP<=pCut]
		nodesP<--log10(nodesP)
		temp1<-table(resultUsed[,1])
		temp1<-temp1[names(nodesP)]
		temp1[is.na(temp1)]<-0
		temp2<-table(resultUsed[,2])
		temp2<-temp2[names(nodesP)]
		temp2[is.na(temp2)]<-0
		nodesCount<-temp1+temp2
		nodes<-cbind(nodesP,nodesCount)
		network<-resultUsed[which(resultUsed[,1] %in% names(nodesP) & resultUsed[,2] %in% names(nodesP)),c(1,2,4)]
		network[,3]<--log10(network[,3])
		networkList<-list(nodes=nodes,network=network)
	} else { #correlationAnalysis result
		temp<-apply(resultUsed[,4:6],1,max,na.rm=T)
		resultUsed<-resultUsed[which(temp<=pCut),]
		temp1<-sapply(split(resultUsed,as.character(resultUsed[,1]),drop = FALSE),function(x) max(x[,5]))
		temp2<-sapply(split(resultUsed,as.character(resultUsed[,2]),drop = FALSE),function(x) max(x[,6]))
		nodesP<-sapply(split(c(temp1,temp2),names(c(temp1,temp2))),max)
		nodesP<--log10(nodesP)
		temp1<-table(resultUsed[,1])
		temp1<-temp1[names(nodesP)]
		temp1[is.na(temp1)]<-0
		temp2<-table(resultUsed[,2])
		temp2<-temp2[names(nodesP)]
		temp2[is.na(temp2)]<-0
		nodesCount<-temp1+temp2
		nodes<-cbind(nodesP,nodesCount)
		network<-resultUsed[,c(1,2,4)]
		network<-resultUsed[which(resultUsed[,1] %in% names(nodesP) & resultUsed[,2] %in% names(nodesP)),c(1,2,4)]
		network[,3]<--log10(network[,3])
		networkList<-list(nodes=nodes,network=network)
	}

	#networkVis
#	require(igraph)
	g <- graph.data.frame(networkList$network, directed=F, vertices=cbind(row.names(networkList$nodes),networkList$nodes))
#	V(g)$size<-data2range(as.numeric(V(g)$"nodesP"),5,20)
	cateNodeSize<-data2cate(as.numeric(V(g)$"nodesP"),cate=c(5,10,15,20))
	V(g)$size<-cateNodeSize$cate
	E(g)$width<-data2range(E(g)$"p.shareTargets",1,5)
	plot.igraph(g,vertex.frame.color=NA,vertex.label.degree=pi,vertex.label.dist=1)
#	legend("topright",legend=cateNodeSize$levels,bty="n",pch=16,cex=1:4)
	return(networkList)
}

data2range<-function(x,minNew=0,maxNew=1) {
	(x-min(x))/(max(x)-min(x))*(maxNew-minNew)+minNew
}

data2cate<-function(x,breaks,cate) {
	if (missing(breaks)) {
		breaks<-c(min(x),min(x)+(max(x)-min(x))*1/4,min(x)+(max(x)-min(x))*2/4,min(x)+(max(x)-min(x))*3/4,max(x))
	}
	temp<-cut(x,breaks=breaks,include.lowest=T)
	result<-cate[temp]
	return(list(cate=result,levels=levels(temp)))
}

pkgInfo<-function(...) {
	cat(paste0("[tfTargetIdentification] ",...))
	cat("\n")
}

targetCor<-function(corMatrix,targets,n) {
	targetsInd<-which(colnames(corMatrix) %in% targets)
	if (n<length(targetsInd)) {
		temp<-sample(targetsInd,n)
	} else {
		temp<-targetsInd
	}
	temp<-unlist(sapply(1:(length(temp)-1),function(x) (temp[x]-1)*ncol(corMatrix)+temp[(x+1):length(temp)]))
	result<-corMatrix[temp]
	return(result)
}