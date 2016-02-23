##' DefineStartAndStop: Define start and stop for a given CNV. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineStartAndStop
##' @return Data frame with start and stop.
##' @author Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples

DefineStartAndStop <- function(indx, subCNV, CHR, ID, CPT.Res)
{
	suppressPackageStartupMessages(library(pastecs))
	
	Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
	Info <- t(Info)
	DF <- as.data.frame(Info)
	CNVMean <- apply(DF, 1, function(X){ mean(CPT.Res@data.set[X[1]:X[2]]) } )
	DF$CNVMean <- CNVMean
	names(DF) <- c("Start", "Stop", "CNVMean")
	# removed many of the controls (normalization should take care of it).

	# Adding index SNP chip position.
	if(sum(DF$Start < 1) > 0){ DF$Start[DF$Start < 1] <- 1 }
	if(sum(DF$Stop > nrow(subCNV)) > 0){ DF$Stop[DF$Stop > nrow(subCNV)] <- nrow(subCNV)}

	DF$StartIndx <- DF$Start
	DF$StopIndx <- DF$Stop
	DF$Start <- subCNV$Position[DF$Start]
	DF$Stop <- subCNV$Position[DF$Stop]
	DF$Length <- DF$Stop - DF$Start
	DF$NumSNPs <- DF$StopIndx - DF$StartIndx
	DF$Chr <- rep(unique(subCNV$Chr), nrow(DF))
	DF$ID <- rep(ID, nrow(DF))
	return(DF)
}
