##' FindCNV.V4.: Find CNVs using Log R Ratio (LRR) data by cpt.meanvar method. 
##'
##' @title FindCNV.V4.
##' @param ID: Sample ID. 
##' @param Sample: A data frame with samples information like LRR and BAF.
##' @param CPTmethod: Unknown, default = meanvar.
##' @param CNVSignal: Unknown, default = 0.1.
##' @param Penvalue: Unknown, defualt = 20.
##' @param Merge: Unknown, default = TRUE.
##' @return A data frame with CNVs, with start and stop position, chromosome (chr.), and etc. for each sample. 
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @export
##' @examples Unknown.

FindCNV.V4 <- function(ID="Test", Sample=Sample, CPTmethod="HMM", CNVSignal=0.1, penvalue=16, Merge=TRUE, minseglen=20)
{
	
	suppressPackageStartupMessages(library("changepoint"))
	suppressPackageStartupMessages(library("depmixS4"))

	tmp <- sapply(unique(Sample$Chr), function(X) # X is chr. Loop over Chr.
	{
		CHR <- as.character(X)
		subSample <- subset(Sample, Chr %in% CHR)
		subSample <- subSample[with(subSample, order(Position)),]
		
		# Getting prob for each SNP using depmix HMM and viterbi. 
		LRR.mod <- depmix(Log.R.Ratio ~ 1, family = gaussian(), nstates = 3, data = subSample, instart=c(0.1, 0.8, 0.1), respstart=c(-0.45,0,0.3, 0.2,0.2,0.2))
		MyFit2 <- setpars(LRR.mod, getpars(HMM.LRR.fit))
		LRR.probs <- viterbi(MyFit2)
		
		# Using changepoint package	
		if(CPTmethod %in% "HMM")
		{
			State <- LRR.probs$state
			indx <- sapply(1:(length(State)-1), function(i){ if(State[i] != State[(i+1)]){ return(i) }})
			indx <- unlist(indx)
		}	
		else if(CPTmethod %in% "meanvar")
		{
			CPT.Res <- cpt.meanvar(subSample$Log.R.Ratio, method='PELT', class=TRUE, pen.value=penvalue, penalty="Manual", minseglen=minseglen)
			indx <- cpts(CPT.Res)
		}
		else if(CPTmethod %in% "mean")
		{
			CPT.Res <- cpt.mean(subSample$Log.R.Ratio, method='PELT', penalty="AIC", class=TRUE, minseglen=minseglen)
			indx <- cpts(CPT.Res)
		}
		
		indx <- c(1, indx, length(subSample$Log.R.Ratio))

		DF <- DefineStartAndStop(indx, subSample, CHR, ID)

		# Adding mean probability of each state for each CNV
		Probs <- apply(DF, 1, function(X){ res <- apply(LRR.probs[as.numeric(X["StartIndx"]):as.numeric(X["StopIndx"]), 2:4], 2, mean); sort(res, decreasing=TRUE)[1] }) 
		DF$prob <- Probs
					 
		# Giving a preliminary CN for each segment
		DF$CN <- DF$CNVMean
		DF$CN[DF$AAB > DF$Het || DF$ABB > DF$Het & DF$CN > 0] <- 3
		DF$CN[DF$CN < 0 & DF$Het < 5] <- 1
		DF$CN[DF$CN != 1 & DF$CN != 3] <- 2
		DF <- subset(DF, abs(CNVMean) > CNVSignal)
		
		if(nrow(DF) > 1) # Changed the pen.value for cpt.meanvar and it does not break much. Maybe no need mergeing.
		{
			if(Merge)
			{
				DF2 <- MergeCNVs(DF, MaxNumSNPs=50)
				DF <- DF2[,colnames(DF)] # returning to same colname order as if you do not go in MergeCNVs
			}
		}
		return(DF)
	})

	df <- MatrixOrList2df(tmp)
	if(nrow(df) > 0)
	{
		df <- df[,!colnames(df) %in% ".id"]
		return(df)
	}
}
