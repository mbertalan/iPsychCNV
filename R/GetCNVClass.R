##' GetCNVClass: read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetCNVClass
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

GetCNVClass <- function(Variable) # Can be 0, 1, 3, 4
{
	Variable <- as.numeric(Variable)
	if(Variable == 0){ Class = "DoubleDel" }
	if(Variable == 1){ Class = "Del" }
	if(Variable == 3){ Class = "Dup" }
	if(Variable == 4){ Class = "DoubleDup" }
	return(Class)
}	
