##' NormalizeData: Normalize LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title Normalize Data
##' @return LRR normalized
##' @author Marcelo Bertalan
##' @export
NormalizeData <- function(CNV,ExpectedMean=0, DF=NA, NormQspline=FALSE)
{
	tmp <- sapply(unique(CNV$Chr), function(X) # X is chr. Loop over Chr.
	{
		#cat(X,"\n")
		subCNV <- CNV[CNV$Chr %in% X,]
		subCNV <- subCNV[with(subCNV, order(Position)),]

		iNorm <- subCNV$Log.R.Ratio >= -1	
		iUnder <- subCNV$Log.R.Ratio < -1
		LRR_SD <- sd(subCNV$Log.R.Ratio)

		if(sum(iUnder) > 5)
		{
			Norm <- subCNV$Log.R.Ratio
			if(LRR_SD > 0.25)
			{
				Norm[iNorm] <-  (subCNV$Log.R.Ratio[iNorm] + (ExpectedMean - mean(subCNV$Log.R.Ratio[iNorm])))/(sd(subCNV$Log.R.Ratio[iNorm])/0.25)
				Norm[iUnder] <-  (subCNV$Log.R.Ratio[iUnder] + (ExpectedMean - mean(subCNV$Log.R.Ratio[iUnder])))/(sd(subCNV$Log.R.Ratio[iUnder])/0.25)
			}
			else
			{
				Norm[iNorm] <-  (subCNV$Log.R.Ratio[iNorm] + (ExpectedMean - mean(subCNV$Log.R.Ratio[iNorm])))
				Norm[iUnder] <-  (subCNV$Log.R.Ratio[iUnder] + (ExpectedMean - mean(subCNV$Log.R.Ratio[iUnder])))
			}
		}
		else
		{
			if(LRR_SD > 0.25)
			{
				Norm <-  (subCNV$Log.R.Ratio + (ExpectedMean - mean(subCNV$Log.R.Ratio)))/(sd(subCNV$Log.R.Ratio)/0.25)
			}
			else
			{
				Norm <-  (subCNV$Log.R.Ratio + (ExpectedMean - mean(subCNV$Log.R.Ratio)))
			}		
		}
		
		subCNV$Log.R.Ratio <- Norm	
		subCNV$Log.R.Ratio <- NormalizeLRR(subCNV$Log.R.Ratio, ExpectedMean=0, DF, NormQspline)
		return(subCNV)
	})
	tmp2 <- MatrixOrList2df(tmp)
	return(tmp2)
}
		

		
