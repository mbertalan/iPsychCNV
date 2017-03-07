##' FilterFromCNVs: Filter CNV from other methods.
##'
##' Specifically designed to reduce false positive CNVs and handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title FilterFromCNVs
##' @param CNVs: Data frame with CNVs.
##' @param PathRawData: The path to the raw data files containing Log R Ratio (LRR) and B Allele Frequency (BAF) values. Example: "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/Version2".
##' @param MinNumSNPs: Minimum number of SNPs per CNV, default = 10. 
##' @param Source: Which method is the original call, default = iPsychCNV.
##' @param Skip: Integer, the number of rows of the data file that should be skipped before beginning tp read the data. Use if file has comments, default = 0.
##' @param Cores: Number of cores to run in parallel, default = 1. 
##' @return Data frame with the estimated copy number for each chromosome.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples 
##' LongRoi <- MakeLongMockSample(Mean=c(-0.6, -0.3, 0.3, 0.6), Size=c(200, 400, 600))
##' # GADA
##' Sample <- read.table("LongMockSample.tab", sep="\t", header=TRUE, stringsAsFactors=F)
##' Gada <- RunGada(Sample)
##' Gada.filter <- FilterFromCNVs(CNVs=Gada, PathRawData=".", MinNumSNPs=10, Source="Gada", Skip=0, Cores=1)
##' See iPsychCNV tutorial for more examples
##' http://biopsych.dk/iPsychCNV/tutorial.html

FilterFromCNVs <- function(CNVs, PathRawData, MinNumSNPs=10, Skip=0, Cores=1, IndxPos=FALSE)
{
	suppressPackageStartupMessages(library(parallel))
	suppressPackageStartupMessages(library("depmixS4"))
	
	if(IndxPos)
	{
		CNVs <- GetIndxPositionFromChips(CNVs, Sample)
	}
	
	tmp <- mclapply(unique(CNVs$ID), mc.cores=Cores, mc.preschedule = FALSE, function(IDs) 
	{
		subCNVs <- subset(CNVs, ID %in% IDs)
		RawFile <- paste(PathRawData, "/", IDs, sep="", collapse="")
		cat(RawFile, "\r")
		Sample <- ReadSample(RawFile, skip=Skip, LCR=FALSE)
		tmp2 <- FilterCNVs.V5(CNVs=subCNVs, MinNumSNPs=MinNumSNPs, Sample=Sample, ID=IDs)
		return(tmp2)
	})
	
	tmp2 <- MatrixOrList2df(tmp)
	df <- tmp2[, !colnames(tmp2) %in% ".id"]
	return(df)
}
