##' BindMultipleObjects: Read files and combine. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title BindMultipleFiles
##' @param PathRawData: Full path where the files are located.
##' @param Pattern: The patter of file name. For example, *.tab
##' @param recursive: If should search for sub-folders. 
##' @return Return data in a data frame.
##' @author Marcelo Bertalan.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples 
##' All.CNV.Files <- BindMultipleFiles(PathRawData=".", Pattern="CNV.tab", recursive=TRUE)

BindMultipleFiles <- function(PathRawData=".", Pattern="CNV.tab", recursive=TRUE)
{
	library(data.table)
	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
	
	tmp2 <- sapply(Files, function(X)
	{
		fread(X, sep="\t", header=TRUE, stringsAsFactor=F) 
	})

	tmp3 <- MatrixOrList2df(tmp2)
	return(tmp3)
}
