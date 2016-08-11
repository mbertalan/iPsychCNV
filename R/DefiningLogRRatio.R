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
	if(res2$CNVmeanByRef < -0.2)
	{
		LogRRatio <- 1
	}
	else if(res2$CNVmeanByRef > 0.05)
	{
		LogRRatio <- 3
	}
	return(LogRRatio)
}
#	if(res2$SDCNV < 0.7)
#	{
#		if(res2$CNVmeanByRef < 0) # 1 or 0
#		{
#			if(res2$CNVmean < res2$HighMean &  res2$CNVmean < res2$LowMean)
#			{
#				if(res2$DiffHigh > 0.7 || res2$DiffLow > 0.7)
#				{
#					if(res2$CNVmean < -1){ LogRRatio <- 0 }
#				}
#				else if(res2$CNVmeanByRef < -0.2)
#				{
#					LogRRatio <- 1
#				}
#				else
#				{
#					LogRRatio <- 1
#				}
#			}
#		}
#		else if(res2$CNVmeanByRef > 0)
#		{
#			if(res2$CNVmean > res2$HighMean &  res2$CNVmean > res2$LowMean)
#			{
#				if(CNVmeanByRef > 0.4)
#				{
#					LogRRatio <- 4
#				}
#				else
#				{
#					LogRRatio <- 3
#				}
#			}
#		}	
#	}		
#	return(LogRRatio)
#}
