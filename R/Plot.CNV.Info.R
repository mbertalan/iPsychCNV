##' Plot.CNV.Info: Define Start and Stop for CNV. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title Plot.CNV.Info
##' @return Data frame with Start and Stop.
##' @author Marcelo Bertalan
##' @export

Plot.CNV.Info <- function(MinNumSNPs, DF, subCNV, ID)
{
	if(sum(DF$Start < 1) > 0){ DF$Start[DF$Start < 1] <- 1 }
	if(sum(DF$Stop > nrow(subCNV)) > 0){ DF$Stop[DF$Stop > nrow(subCNV)] <- nrow(subCNV)}

	DF$StartIndx <- DF$Start
	DF$StopIndx <- DF$Stop
	DF$Start <- subCNV$Position[DF$Start]
	DF$Stop <- subCNV$Position[DF$Stop]
	DF$Length <- DF$Stop - DF$Start
	DF$NumSNPs <- DF$StopIndx - DF$StartIndx
	DF <- subset(DF, NumSNPs > MinNumSNPs)
	DF$Chr <- rep(unique(subCNV$Chr), nrow(DF))
	DF$ID <- rep(ID, nrow(DF))
	return(DF)
}
