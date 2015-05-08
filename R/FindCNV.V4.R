##' FindCNV: Find CNVs using Log.R.Ratio data by cpt.meanvar method. 
##'
##' Input: sample path, CNV in data frame format and min number of SNPs
##' @title FindCNV
##' @return A data frame with CNVs, with start and stop position, chr and etc. for each sample. 
##' @author Marcelo Bertalan
##' @export

FindCNV.V4 <- function(ID, MinNumSNPs, CNV, CPTmethod="meanvar", CNVSignal=0.19)
{
	
	suppressPackageStartupMessages(library("changepoint"))

	tmp <- sapply(unique(CNV$Chr), function(X) # X is chr. Loop over Chr.
	{
		CHR <- as.character(X)
		subCNV <- subset(CNV, Chr %in% CHR)
		subCNV <- subCNV[with(subCNV, order(Position)),]
		
		# Using changepoint package	
		if(CPTmethod %in% "meanvar")
		{
			CPT.Res <- cpt.meanvar(subCNV$Log.R.Ratio, method='PELT', class=TRUE)
		}
		else if(CPTmethod %in% "mean")
		{
			CPT.Res <- cpt.mean(subCNV$Log.R.Ratio, method='PELT', penalty="AIC", class=TRUE)
		}
		indx <- cpts(CPT.Res)
		indx <- c(1, indx, length(subCNV$Log.R.Ratio))

		DF <- DefineStartAndStop(indx, subCNV, MinNumSNPs, CHR, ID, CPT.Res)
		
		# Using meanvar it breaks CNVs, I am trying to merge it again.
		DF <- subset(DF, abs(CNVMean) > 0.1)
		DF$CN <- DF$CNVMean
		DF$CN[DF$CN > 0] <- 3
		DF$CN[DF$CN < 0] <- 1
		DF2 <- MergeCNVs(DF)
		
		#Df <- FixCNVPosition(DF, subCNV, MinNumSNPs, ID) # Fix the position of the CNV by distance of the SNPs
		return(DF2)
	})

	df <- MatrixOrList2df(tmp)
	if(nrow(df) > 0)
	{
		df <- subset(df, NumSNPs > MinNumSNPs)
		df <- df[,!colnames(df) %in% ".id"]
		return(df)
	}
}

