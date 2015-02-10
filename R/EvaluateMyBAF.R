##' EvaluateMyBAF: Evaluate BAF. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title EvaluateMyBAF
##' @return Classification for BAF.
##' @author Marcelo Bertalan
##' @export

EvaluateMyBAF <- function(res, CNVmean)
{
	#Centers <- c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1)
	#names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	
	if(CNVmean < -0.15)
	{
		if(res$AAAA > 15  & res$AB < 8 & res$BBBB > 15 & res$ABB < 9 & res$AAB < 6)  # Del
		{
			BAlleleFreq <- "1"
		}
		else if(res$AAAA > 8 & res$BBBB > 8 & res$AB > 8 & res$AAB > 8 & res$ABB > 8 & res$AAAB > 8 & res$ABBB > 8) # Double del
		{
			BAlleleFreq <- "0"
		}
		else
		{
			BAlleleFreq <- "Undefined"
		}
	}
	else if(CNVmean > 0.15)
	{
		if(res$AAB > 3 & res$AB < 6 & res$ABB > 3) # Dup
		{
			BAlleleFreq <- "3"
		}
		else if(res$AAAB > 4 & res$ABBB > 4 & res$AB > 4) # Double Dup AB = AABB
		{
			BAlleleFreq <- "4"
		}
		else
		{
			BAlleleFreq <- "Undefined"
		}
	}
	else
	{
		BAlleleFreq <- "Undefined"		
	}
	return(BAlleleFreq)
}
