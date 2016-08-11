##' MatrixOrList2df: Read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MatrixOrList2df
##' @param tmp: Unknwon, default = Unknown.
##' @return return data in a data frame.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

MatrixOrList2df <- function(tmp)
{
	library(plyr)
	library(data.table)
	if(class(tmp) %in% "matrix")
	{
		df <- apply(tmp, 2, function(X){ as.data.frame(X,  stringsAsFactors=F) })
		#df <- ldply (df, data.frame)
		df <- do.call(rbind.data.frame, df)
	}
	if(class(tmp) %in% "list")
	{	
		#tmp <- tmp[unlist(lapply(tmp, function(X){ names(X)[1] %in% "Start" }))]
		#df <- do.call(rbind.data.frame, tmp)
		Indx <- lapply(tmp, is.data.frame)
		tmp2 <- tmp[unlist(Indx)]
		tmp3 <- lapply(tmp2, colnames) # Get colnames from all elements in the list 
		tmp3 <- Reduce(intersect, tmp3) # Get only common colnames. 
		tmp2 <- lapply(tmp2, function(X){ X[,tmp3]}) # Making the list with same columns. 
		df <- rbindlist(tmp2)
	}
	if(class(tmp) %in% "data.frame")
	{
		df <- tmp
	}
	return(df)
}
