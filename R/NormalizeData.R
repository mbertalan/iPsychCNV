##' NormalizeData: Normalize LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title Normalize Data
##' @return LRR normalized
##' @author Marcelo Bertalan
##' @export

NormalizeData <- function(CNV,ExpectedMean=0, penalty=20, Quantile=TRUE, QSpline=TRUE, sd=0.2)
{
	library(preprocessCore)
	tmp <- sapply(unique(CNV$Chr), function(X) # X is chr. Loop over Chr.
	{
		subCNV <- CNV[CNV$Chr %in% X,]
		subCNV <- subCNV[with(subCNV, order(Position)),]
		subCNV$LRR <- subCNV$Log.R.Ratio
		LRR <- subCNV$Log.R.Ratio
		
		# add later data that LRR with 2 peaks.
		
		if(QSpline) # detrend the data, only when sd is high & sd(LRR) > 0.2
		{
			Spline <- smooth.spline(LRR, penalty=penalty)
			Mean <- Spline$y
			LRR <- LRR - Mean
		}
	
		if(Quantile) # Same distribution, fixed sd and mean
		{
			# Creating perfect data
			M <- sapply(1:50, function(N){ rnorm(n=length(LRR), mean=0, sd=sd) })
			M2 <- cbind(M, LRR)
			M3 <- normalize.quantiles(M2)
			LRR <- M3[, 51]
		}
		subCNV$Log.R.Ratio <- LRR
		return(subCNV)
	})
	tmp2 <- MatrixOrList2df(tmp)
	return(tmp2)
}
