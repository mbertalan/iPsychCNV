##' MatrixOrList2df: read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MatrixOrList2df
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

MatrixOrList2df <- function(tmp)
{
	library(plyr)
	if(class(tmp) %in% "matrix")
	{
		df <- apply(tmp, 2, function(X){ as.data.frame(X,  stringsAsFactors=F) })
		#df <- ldply (df, data.frame)
		df <- do.call(rbind.data.frame, df)
	}
	if(class(tmp) %in% "list")
	{	
		#tmp <- tmp[unlist(lapply(tmp, function(X){ names(X)[1] %in% "Start" }))]
		df <- do.call(rbind.data.frame, tmp)
	}
	if(class(tmp) %in% "data.frame")
	{
		df <- tmp
	}
	return(df)
}
