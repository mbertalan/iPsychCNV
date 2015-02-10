##' FixCNVPosition: Trim the CNV position by distance of SNPs. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title FixCNVPosition
##' @return Data frame with CNV information.
##' @author Marcelo Bertalan
##' @export

FixCNVPosition <- function(Df, subCNV, MinNumSNPs, ID)
{
	tmp <- apply(Df, 1, function(X)
	{
		Start <- as.numeric(X["StartIndx"])
		Stop <- as.numeric(X["StopIndx"])
		SNPsDistance <- subCNV$Position[(Start+1):(Stop)] - subCNV$Position[(Start):(Stop-1)]
		SNPsDistance[6:(length(SNPsDistance)-6)] <- 100 # Avoid to break the CNV in half. Just want to cut the edges. 
		Newindx <- SNPsDistance > 100000 # mean(SNPsDistance)*4
		Newindx <- which(Newindx)
		Newindx <- Newindx + (Start - 1)
		Newindx <- c(Start, Newindx, Stop)
		return(Newindx)
	})
	tmp2 <- unlist(tmp)
	Info <- sapply(1:(length(tmp2)-1), function(X){ rbind(tmp2[X],tmp2[(X+1)]) })
	Info <- t(Info)
	Df <- as.data.frame(Info)
	colnames(Df) <- c("Start", "Stop")
	DF <- Plot.CNV.Info(MinNumSNPs, Df, subCNV, ID)
	return(DF)
}
