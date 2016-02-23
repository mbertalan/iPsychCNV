##' DefineBAFType: Define B Allele Frequency (BAF) by peaks. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineBAFType
##' @param SumPeaks: Unknown
##' @return BAF classification.
##' @author Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown

DefineBAFType <- function(SumPeaks)
{
	if(SumPeaks == 4 ) # | ClusterNum == 4
	{
		BAlleleFreq <- 3
	}
	else if(SumPeaks == 5) # | ClusterNum == 4
	{
		BAlleleFreq <- 4
	}
	else if(SumPeaks == 2 ) # | ClusterNum == 2
	{
		BAlleleFreq <- 1
	}
	else
	{
		BAlleleFreq <- "Undefined"
	}
	df <- data.frame(BAlleleFreq=BAlleleFreq, stringsAsFactors=FALSE)
	return(df)
}
