##' DefineCNVType: Define if the CNV is good or bad. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineCNVType
##' @return CNV type.
##' @author Marcelo Bertalan
##' @export

DefineCNVType <- function(tmp2)
{
	Type <- apply(tmp2, 1, function(X)
	{
		if(X["LogRRatio"] %in% X["BAlleleFreq"] & X["BAlleleFreq"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "Undefined" )
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
