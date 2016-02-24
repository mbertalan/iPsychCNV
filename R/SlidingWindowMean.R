##' SlideWindowMean: Create a smooth effect on the data. 
##'
##' X <- a vector with the data. N <- the window size. 
##' @title SlideWindowMean
##' @param X: Unknown, default = Unknown.
##' @param Y: Unknown, default = Unknown. 
##' @return A vector with mean of the window size.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown. 
##' 

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
