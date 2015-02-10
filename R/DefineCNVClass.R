##' DefineCNVClass: Define if the CNV class. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineCNVClass
##' @return CNV class.
##' @author Marcelo Bertalan
##' @export


DefineCNVClass <- function(tmp2)
{
	Type <- apply(tmp2, 1, function(X)
	{
		Class <- "Normal"
		if(X["LogRRatio"] %in% X["BAlleleFreq"] & X["BAlleleFreq"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "Undefined" )
		{
			Class <- GetCNVClass(X["LogRRatio"])
		}
		else if(X["LogRRatio"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "Undefined" ) # Just to avoid 2 Undefined being good.
		{	
				
			Class <- GetCNVClass(X["LogRRatio"])
		}
		else if(X["LogRRatio"] %in% X["BAlleleFreq"] & !X["LogRRatio"] %in% "Undefined")
		{	
			Class <- GetCNVClass(X["LogRRatio"])
		}
		return(Class)
	})
	return(Type)
}
