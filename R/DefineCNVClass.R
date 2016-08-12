##' DefineCNVClass: Define the CNV class. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineCNVClass
##' @param tmp2: Unknown
##' @return CNV class.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Uknown
##'

DefineCNVClass <- function(tmp2)
{
	Type <- apply(tmp2, 1, function(X)
	{
		Class <- 2
		if(X["LogRRatio"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "2" )
		{
			#Class <- GetCNVClass(X["LogRRatio"])
			Class <- X["LogRRatio"]
		}
		else if(X["LogRRatio"] %in% X["MyBAF"] & !X["LogRRatio"] %in% "2" ) # Just to avoid 2 Undefined being good.
		{	
				
			#Class <- GetCNVClass(X["LogRRatio"])
			Class <- X["LogRRatio"]
		}
		else
		{
			if(X["LogRRatio"] == 1 & X["MyBAF"] == 0)
			{
				Class <- 0
			}
			else if(X["LogRRatio"] == 3 & X["MyBAF"] == 4)
			{
				Class <- 4
			}
			else if(X["LogRRatio"] == 3 & X["MyBAF"] == 1) # Duplication with LOH.
			{
				Class <- 2 # removing 5, not working well.
			}
		}
		
		return(Class)
	})
	return(Type)
}
