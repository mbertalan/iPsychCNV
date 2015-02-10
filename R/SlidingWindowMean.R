##' SlideWindowMean: Create a smooth effect on the data. 
##'
##' X <- a vector with the data. N <- the window size. 
##' @title SlideWindowMean
##' @return A vector with mean of the window size.
##' @author Marcelo Bertalan
##' @export

SlideWindowMean <- function(X, N)
{
	#N <- round(length(X)/10)
	sapply(1:length(X), function(Z)
	{
		if((Z+N) <= length(X))
		{
			mean(X[Z:(Z+N)]) 
		}
		else
		{
			mean(X[Z:length(X)])
		}
	})
}
