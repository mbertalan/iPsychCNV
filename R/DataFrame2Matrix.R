##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##' Specifically designed to reduce false positive CNVs and amplified DNA on dried blood spots.

##' DataFrame2Matrix: Transform 3 vectors from a data frame into matrix.
##'
##' @title DataFrame2Matrix

##' @param RowNames Vector that will be matrix's rownames.
##' @param ColNames Vector that will be matrix's colnames.
##' @param Values Vector that will be matrix's values.
##' @return return a matrix with values and specific rownames and colnames.

##' @author Marcelo Bertalan
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' M <- DataFrame2Matrix(RowNames=V1, ColNames=V2, Values=V3)

DataFrame2Matrix <- function(RowNames=V1, ColNames=V2, Values=V3)
{
	DF <- data.frame(V1=RowNames, V2=ColNames, V3=Values)

	M <- matrix(nrow=length(unique(RowNames)), ncol=length(unique(ColNames)))
	colnames(M) <- unique(ColNames)
	rownames(M) <- unique(RowNames)

	apply(DF, 1, function(X)
	{
		M[X[1],X[2]] <<- as.numeric(X[3])
	})

	return(M)
}
