##' DefineStartAndStop: Define Start and Stop for CNV. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineStartAndStop
##' @return Data frame with Start and Stop.
##' @author Marcelo Bertalan
##' @export

DefineStartAndStop <- function(indx, subCNV, MinNumSNPs, CHR, ID, CPT.Res)
{
	suppressPackageStartupMessages(library(pastecs))
	
	Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
	Info <- t(Info)
	Df <- as.data.frame(Info)
	CNVMean <- apply(Df, 1, function(X){ mean(CPT.Res@data.set[X[1]:X[2]]) } )
	Df$CNVMean <- CNVMean
	names(Df) <- c("Start", "Stop", "CNVMean")
	# removed many of the controls (normalization should take care of it).

	DF3 <- Plot.CNV.Info(MinNumSNPs, Df, subCNV, ID) # Min number of SNPs for CNV
	
	return(DF3)
}
