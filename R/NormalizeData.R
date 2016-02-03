##' NormalizeData: Normalize LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title Normalize Data
##' @return LRR normalized
##' @author Marcelo Bertalan
##' @export

NormalizeData <- function(Sample=Sample,ExpectedMean=0, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18)
{
	library(preprocessCore)
	tmp <- sapply(unique(Sample$Chr), function(X) # X is chr. Loop over Chr.
	{
		subSample <- Sample[Sample$Chr %in% X,]
		subSample <- subSample[with(subSample, order(Position)),]
		subSample$LRR <- subSample$Log.R.Ratio
		LRR <- subSample$Log.R.Ratio
	
		# Center chromosome LRR	
		subSample$Log.R.Ratio <- subSample$Log.R.Ratio - mean(subSample$Log.R.Ratio)
		
		# add later data that LRR with 2 peaks.
		if(QSpline) # detrend the data, only when sd is high & sd(LRR) > 0.2
		{
			Spline <- smooth.spline(LRR, penalty=penalty)
			Mean <- Spline$y
			LRR <- LRR - Mean
			subSample$Log.R.Ratio <- LRR
		}
	
		if(Quantile) # Same distribution, fixed sd and mean
		{
			# Creating perfect data
			M <- sapply(1:50, function(N){ rnorm(n=length(LRR), mean=0, sd=sd) })
			M2 <- cbind(M, LRR)
			M3 <- normalize.quantiles(M2)
			LRR <- M3[, 51]
			subSample$Log.R.Ratio <- LRR
		}
		
		return(subSample)
	})
	tmp2 <- MatrixOrList2df(tmp)
	return(tmp2)
}
