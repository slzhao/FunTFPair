FunTFPair
============
* [Introduction](#Introduction)
* [Required packages](#require)
* [Download and install](#download)
* [Browse vignette](#vignette)
* [Examples](#example)

<a name="Introduction"/>
# Introduction #
Summary: Transcription factors (TFs) are fundamental regulators of gene expression and generally function in a complex and cooperative manner. Identifying context-dependent cooperative TFs is essential for understanding how cells respond to environmental change. The huge amount of various omics data currently available, providing genome-wide physical binding and functional effect information about transcription, holds a great opportunity to study TF cooperativity. Here we developed an R package, FunTFPair, which provides an easy and powerful way to identify condition-specific TF pairs by integrating transcription factor binding sites from ENCODE and gene expression profiles from GEO. Users only need to provide GEO ID that they are interested in. FunTFPair will automatically retrieve the expression profiles of the input GEO ID, get TF targets from ENCODE, and identify TF pairs whose common targets show statistically significant differences under changed experimental conditions or exhibit coordinated transcription in a specific condition. The functional TF pairs and their relative importance will be reported in a TF cooperation network. Two datasets from GEO are used as examples to demonstrate the usage and the reliability of the package.


Availability and Implementation: FunTFPair is implemented in R language.  The source code is provided at https://github.com/slzhao/FunTFPair. A more detailed vignette is also released along the package.

<a name="require"/>
# Required packages #
FunTFPair package requires some other R packages, including limma, biomaRt, GEOquery and igraph. If you haven't installed them in your computer, you can use the following R codes to install them:
	
	#Running the following codes in your R	
	install.packages("igraph")
	source("http://bioconductor.org/biocLite.R")
	biocLite(pkgs=c("limma","biomaRt","GEOquery"))

<a name="download"/>
# Download and install #
You can download and install FunTFPair package from [github](https://github.com/slzhao/FunTFPair) by the following R codes.

	#Running the following codes in your R
	library(devtools)
    install_github("slzhao/FunTFPair")

<a name="vignette"/>
# Browse vignette #
A vignette, which included introduction for functions, executable examples and description for TF cooperative network result,  was also provided with this package. You can use the following R codes to browse it:

	#View the vignette
	browseVignettes("FunTFPair")

<a name="example"/>
# Examples #
After you have installed FunTFPair package. You can enter R and use following R codes to see the examples for it.
	

	#Load package
	library("FunTFPair")
	#Examples for GEO data preparation
	example(prepareGeoData)
	#Examples for differential analysis based functional TF pair identification
	example(differentialAnalysis,run.dontrun=TRUE)
	#Examples for correlation analysis based functional TF pair identification
	example(correlationAnalysis,run.dontrun=TRUE)
	#Examples for functional TF Network visualization
	example(networkVis,run.dontrun=TRUE)
