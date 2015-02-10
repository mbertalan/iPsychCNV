##' AddInfo2res: read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title AddInfo2res
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

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
