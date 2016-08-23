##' Char2Num: Function to transform character columns in data frame to numeric, if all are numeric. 
##' 
##' Specifically designed to reduce false positive CNVs and amplified DNA on dried blood spots.
##' @title Char2Num
##' @param df: Input data frame. 
##' @return Return the data frame with numeric columns if all values in column are numeric.
##' @author Marcelo Bertalan.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' df2 <- Char2Num(df)

Char2Num <- function(df)
{
	Indx <- apply(df, 2, function(X){ sum(is.na(as.numeric(X))) })
	df[,Indx == 0] <- sapply(df[, Indx == 0], as.numeric)
	return(df)
}
