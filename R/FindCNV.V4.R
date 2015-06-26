##' FindCNV: Find CNVs using Log.R.Ratio data by cpt.meanvar method. 
##'
##' Input: sample path, CNV in data frame format and min number of SNPs
##' @title FindCNV
##' @return A data frame with CNVs, with start and stop position, chr and etc. for each sample. 
##' @author Marcelo Bertalan
##' @export

FindCNV.V4 <- function(ID="Test", MinNumSNPs=2, Sample=Sample, CPTmethod="meanvar", CNVSignal=0.1, penvalue=10)
{
	
	suppressPackageStartupMessages(library("changepoint"))

	tmp <- sapply(unique(Sample$Chr), function(X) # X is chr. Loop over Chr.
	{
		CHR <- as.character(X)
		subSample <- subset(Sample, Chr %in% CHR)
		subSample <- subSample[with(subSample, order(Position)),]
		
		# Using changepoint package	
		if(CPTmethod %in% "meanvar")
		{
			CPT.Res <- cpt.meanvar(subSample$Log.R.Ratio, method='PELT', class=TRUE, pen.value=penvalue, penalty="Manual")
		}
		else if(CPTmethod %in% "mean")
		{
			CPT.Res <- cpt.mean(subSample$Log.R.Ratio, method='PELT', penalty="AIC", class=TRUE)
		}
		indx <- cpts(CPT.Res)
		indx <- c(1, indx, length(subSample$Log.R.Ratio))

		DF <- DefineStartAndStop(indx, subSample, MinNumSNPs, CHR, ID, CPT.Res)
		
		# Using meanvar it breaks CNVs, I am trying to merge it again.
		DF <- subset(DF, abs(CNVMean) > 0.14)
		DF$CN <- DF$CNVMean
		DF$CN[DF$CN > 0] <- 3
		DF$CN[DF$CN < 0] <- 1
		
		if(nrow(DF) > 1) # Changed the pen.value for cpt.meanvar and it does not break much. Maybe no need mergeing.
		{
			DF2 <- MergeCNVs(DF, MaxNumSNPs=MinNumSNPs)
			DF <- DF2[,colnames(DF)] # returning to same colname order as if you do not go in MergeCNVs
		}
		return(DF)
	})

	df <- MatrixOrList2df(tmp)
	if(nrow(df) > 0)
	{
		df <- subset(df, NumSNPs > MinNumSNPs)
		df <- df[,!colnames(df) %in% ".id"]
		return(df)
	}
}
