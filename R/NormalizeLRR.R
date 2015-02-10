##' NormalizeLRR: Normalize LRR. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title NormalizeLRR
##' @return LRR normalized
##' @author Marcelo Bertalan
##' @export

NormalizeLRR <- function(X, ExpectedMean=0, DF=40, NormQspline=TRUE) # X is a vector
{
	#cat("DF spline: ", DF, "\n")
	if(NormQspline)
	{
		if(is.na(DF)){ LRR.spline <- smooth.spline(X)$y }
		else{  LRR.spline <- smooth.spline(X, df=DF)$y }
	}
	else
	{
		LRR.spline <- X
	}
	return(LRR.spline)
}
