##' DefiningLogRRatio: Define the Log R Ratio (LRR). 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefiningLogRRatio
##' @param res2: Unknown.
##' @return Classification for LRR.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

DefiningLogRRatio <- function(res2)
{
	LogRRatio <- 2
	if(res2$CNVmeanByRef < -0.2) # res2$DiffHigh > 0.1 &  res2$DiffLow > 0.1
	{
		LogRRatio <- 1
	}
	else if(res2$CNVmeanByRef > 0.05) # res2$DiffHigh > 0.02 &  res2$DiffLow > 0.02
	{
		LogRRatio <- 3
	}
	return(LogRRatio)
}
