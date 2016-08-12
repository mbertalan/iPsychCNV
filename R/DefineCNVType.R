##' DefineCNVType: Determine true CNVs calls. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineCNVType
##' @param tmp2: Unknown. 
##' @return CNV type.
##' @author Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

DefineCNVType <- function(tmp2)
{
	Type <- apply(tmp2, 1, function(X)
	{
		if(X["LogRRatio"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "Undefined" )
		{
			return("Good")
		}
		else if(X["LogRRatio"] %in% X["MyBAF"]  )
		{	
				if(!X["LogRRatio"] %in% "Undefined") # Just to avoid 2 Undefined being good.
				{
					return("Good")
				}
				else
				{
					return("Bad")
				}
		}
		else
		{
			return("Bad")
		}
	})
	return(Type)
}
