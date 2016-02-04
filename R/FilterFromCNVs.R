##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##' Specifically designed to reduce false positive CNVs and handle data from amplified DNA on dried blood spots.
##'
##' FilterFromCNVs: Filter CNV from other methods.
##'
##' @title FilterFromCNVs.
##' @param CNVs data frame with CNVs.
##' @param PathRawData Path where are the LLR and BAF files. Exaple: "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/Version2".
##' @param Cores Number of cores to run in parallel. 
##' @param Skip How many rows should skip. Use if file has comments.
##' @param MinNumSNPs Minimum number of SNPs per CNV. 
##' @param Source Which method is the original call.
##' @return Data frame with the estimate copy number for each chromosome.
##' @author Marcelo Bertalan
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

FilterFromCNVs <- function(CNVs, PathRawData, MinNumSNPs=10, Source="iPsychCNV", Skip=0, Cores=1)
{
	suppressPackageStartupMessages(library(parallel))
	
	tmp <- mclapply(unique(CNVs$ID), mc.cores=Cores, mc.preschedule = FALSE, function(IDs) 
	{
		subCNVs <- subset(CNVs, ID %in% IDs)
		RawFile <- paste(PathRawData, "/", IDs, sep="", collapse="")
		cat(IDs, "\r")
		Sample <- ReadSample(RawFile, skip=Skip, LCR=FALSE)

		Res <- apply(subCNVs, 1, function(X)
		{
			StartM <- as.numeric(X["Start"])
			StopM <- as.numeric(X["Stop"])
			ChrM <- X["Chr"]
			ChrM <- gsub(" ","", ChrM)

			subSample <- subset(Sample, Chr %in% ChrM)
			subSample <- subSample[with(subSample, order(Position)),]
			
			StartIndx <- which(subSample$Position == StartM)[1] 
			StopIndx <- which(subSample$Position == StopM)[1]
			df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx)
			cat(ChrM, StartM, StopM, StartIndx, StopIndx, RawFile, "\r")
			return(df)
		})
		#save(Res, file="Res.RData")
		df <- MatrixOrList2df(tmp=Res)
		tmp <- cbind(subCNVs, df)
		#save(tmp, file="tmp.RData")
		tmp2 <- FilterCNVs.V4(CNVs=tmp, MinNumSNPs=MinNumSNPs, Sample=Sample, ID=IDs)
		save(tmp2, file="tmp2.RData")
		return(tmp2)
	})
	#save(tmp, file="tmp.RData")
	tmp2 <- MatrixOrList2df(tmp)
	df <- tmp2[, !colnames(tmp2) %in% ".id"]
	df$Source <- Source
	return(df)
}
