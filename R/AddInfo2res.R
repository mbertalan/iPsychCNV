##' AddInfo2res: Function to add variables for each CNV region. Function is used at FilterCNVs.V4.R. 
##' 
##' Specifically designed to reduce false positive CNVs and amplified DNA on dried blood spots.
##' @title AddInfo2res
##' @param Res:cbind of two data frames. 1-) GetDataVariables. 2-) ClassNumbers.
##' @param CNV2HighPvalue: P value of t.test among CNV region and its upstream region.
##' @param CNV2LowPvalue: P value of t.test among CNV region and its downstream region.
##' @param Class: Result from DefineBAFType function.
##' @param BAlleleFreq: Result from DefineBAFType function.
##' @param MyBAF: Result from EvaluateMyBAF function.
##' @param LogRRatio: Result from  DefiningLogRRatio function.
##' @param SumPeaks: Result from CleaningPeaks function.
##' @param SDChr: Log R Ratio (LRR) standard deviation from whole chromosome. 
##' @param MeanChr: Log R Ratio (LRR) mean from whole chromosome.
##' @return Return all variables together in a data frame.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' res4 <- AddInfo2res(res3, CNV2HighPvalue, CNV2LowPvalue, Class, BAlleleFreq, MyBAF, LogRRatio, SumPeaks, SDChr, MeanChr)

AddInfo2res <- function(res, CNV2HighPvalue, CNV2LowPvalue, Class, BAlleleFreq, MyBAF, LogRRatio, SumPeaks, SDChr, MeanChr)
{
	res$CNV2HighPvalue <- CNV2HighPvalue	
	res$CNV2LowPvalue <- CNV2LowPvalue
	res$Class <- Class
	res$BAlleleFreq <- BAlleleFreq
	res$LogRRatio <- LogRRatio
	res$SumPeaks <- SumPeaks
	res$MyBAF <- MyBAF
	res$SDChr <- SDChr
	res$MeanChr <- MeanChr
	return(res)
}
