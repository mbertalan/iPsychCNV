##' NormalizeLRR: Normalize LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title NormalizeLRR
##' @return LRR normalized
##' @author Marcelo Bertalan
##' @export

NormalizeLRR <- function(X, ExpectedMean=0, DF=NA, NormQspline=FALSE, Quantile=TRUE) # X is a vector
{
	library(preprocessCore)
	#cat("DF spline: ", DF, "\n")
	if(NormQspline)
	{
		if(is.na(DF)){ LRR <- smooth.spline(X)$y }
		else{  LRR <- smooth.spline(X, df=DF)$y }
	}
	else if(Quantile)
	{
		# Creating perfect data
		M <- sapply(1:50, function(N){ rnorm(n=length(X), mean=0, sd=0.2) })
		M2 <- cbind(M, X)
		M3 <- normalize.quantiles(M2)
		LRR <- M3[, 51]
	else
	{
		LRR <- X
	}
	return(LRR)
}
