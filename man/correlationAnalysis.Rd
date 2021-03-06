% Generated by roxygen2 (4.0.2): do not edit by hand
\name{correlationAnalysis}
\alias{correlationAnalysis}
\title{correlationAnalysis}
\usage{
correlationAnalysis(expression, pairs2Target = pairs2TargetRemDup, n = 300)
}
\arguments{
\item{expression}{the expression matrix, can be the result of \code{\link{prepareGeoData}}.}

\item{pairs2Target}{the targets of tf pairs}

\item{n}{the number of genes selected for prior comparison, to make the analysis faster}
}
\description{
A function to perform correlation analysis for expression data, to see if the shared targets of two tf have significantly higher correlations.
}
\details{
A function to perform correlation analysis for expression data, to see if the shared targets of two tf have significantly higher correlations.
}
\examples{
\dontrun{
#Download GEO dataset
dataMatrix3<-prepareGeoData(GEO="GSE54698",geneSymColumnName="GENE_SYMBOL")
dataMatrix3Used1<-dataMatrix3[,1:12]
dataMatrix3Used2<-dataMatrix3[,c(13:24)]
#Choose the TF pairs in HELA cells
tfPairsUsed<-pairs2TargetRemDup[which(pairs2TargetRemDup[,3]=="HELAS3"),]
#Correlation analysis
correlationResult3New1<-correlationAnalysis(dataMatrix3Used1,pairs2Target=tfPairsUsed,n=100)
correlationResult3New2<-correlationAnalysis(dataMatrix3Used2,pairs2Target=tfPairsUsed,n=100)
#Network visualization for control samples
correlationResultNew1Network<-networkVis(correlationResult3New1)
#Network visualization for infected samples
correlationResultNew2Network<-networkVis(correlationResult3New2)
}
}

