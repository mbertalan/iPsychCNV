##' DefineStartAndStop: Define start and stop positions for a given CNV. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineStartAndStop
##' @param Indx: Unknown.
##' @param SubCNV: Unknown.
##' @param CHR: Character, select a specific chromosome to be analyzed. 
##' @param ID: Unknown.
##' @param CPT.Res: Unknown.
##' @return Data frame with start and stop positions for each CNV.
##' @author Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

DefineStartAndStop <- function(indx, subCNV, CHR, ID)
{
	suppressPackageStartupMessages(library(pastecs))
	
	Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
	Info <- t(Info)
	DF <- as.data.frame(Info)
	# Find mean LRR from the CNV region
	CNVMean <- apply(DF, 1, function(X){ mean(subCNV$LRR[as.numeric(X[1]):as.numeric(X[2])]) } )
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

	# Adding %AB
	B <- subCNV$B.Allele.Freq
	AB <- apply(DF[,c("StartIndx", "StopIndx")], 1, function(X){ (sum(B[X[1]:X[2]] > 0.4 & B[X[1]:X[2]] < 0.6)/length(B[X[1]:X[2]]))*100 })
	DF$Het <- AB
	
	#ABB <- apply(DF[,c("StartIndx", "StopIndx")], 1, function(X){ (sum(B[X[1]:X[2]] > 0.6 & B[X[1]:X[2]] < 0.8)/length(B[X[1]:X[2]]))*100 })
	#AAB <- apply(DF[,c("StartIndx", "StopIndx")], 1, function(X){ (sum(B[X[1]:X[2]] > 0.2 & B[X[1]:X[2]] < 0.4)/length(B[X[1]:X[2]]))*100 })
	#DF$ABB <- ABB
	#DF$AAB <- AAB
	
	# Giving a preliminary CN for each segment
	DF$CN <- DF$CNVMean
	DF$CN[DF$CN > 0] <- 3
	DF$CN[DF$CN < 0 & DF$Het < 5] <- 1
	DF$CN[DF$CN <= 0] <- 2
	return(DF)
}
