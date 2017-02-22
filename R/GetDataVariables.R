##' GetDataVariables: Unknown.
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetDataVariables
##' @param Data: Unknown.
##' @return Hotspots - Unknown.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

GetDataVariables <- function(Data)
{
	CNVmean <- mean(Data$CNV) 
	HighMean <- mean(Data$High)
	LowMean <- mean(Data$Low)
	SDCNV <- sd(Data$LRR)
	DiffHigh <- abs(CNVmean - HighMean)
	DiffLow <- abs(CNVmean - LowMean)
	
	res2 <- data.frame(CNVmean=CNVmean, SDCNV=SDCNV, DiffHigh=DiffHigh, DiffLow=DiffLow)
	return(res2)		
}
