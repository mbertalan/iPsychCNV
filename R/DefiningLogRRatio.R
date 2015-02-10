##' DefiningLogRRatio: Define LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefiningLogRRatio
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

DefiningLogRRatio <- function(res2)
{
	#cat("I am in !\n")
	#cat(res2$SDCNV, res2$CNVmean, "\n")
	LogRRatio <- "Undefined"
	if(res2$SDCNV < 0.7)
	{
		if(res2$CNVmean < -0.15) # 1 or 0
		{
			LogRRatio <- "1"
			if(res2$CNVmean < res2$HighMean &  res2$CNVmean < res2$LowMean)
			{
				LogRRatio <- "1"
				if(res2$DiffHigh > 0.7 || res2$DiffLow > 0.7)
				{
					if(res2$CNVmean < -2){ LogRRatio <- "0" }else{ LogRRatio <- "1" }
				}
				else
				{
					LogRRatio <- "1"
				}
			}
		}
		else if(res2$CNVmean > 0.15)
		{
			LogRRatio <- "3"
			if(res2$CNVmean > res2$HighMean &  res2$CNVmean > res2$LowMean)
			{
				LogRRatio <- "3"
				if(res2$DiffHigh > 0.8 || res2$DiffLow > 0.8 & res2$CNVmean > 2)
				{
					LogRRatio <- "4"
				}
				else
				{
					LogRRatio <- "3"
				}
			}
		}	
	}		
	return(LogRRatio)
}
