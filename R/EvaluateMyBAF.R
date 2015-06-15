##' EvaluateMyBAF: Evaluate BAF. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title EvaluateMyBAF
##' @return Classification for BAF.
##' @author Marcelo Bertalan
##' @export

EvaluateMyBAF <- function(res, res2)
{
	#Centers <- c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1)
	#names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	
	if(res2$CNVmeanOrig < -0.1) # It might be the whole chromosome a CNV
	{
		if(res$AAAA > 8  & res$AB < 5 & res$BBBB > 8 & res$ABB < 9 & res$AAB < 9)  # Del
		{
			BAlleleFreq <- 1
		}
		else if(res$AAAA > 5 & res$BBBB > 5 & res$AB > 8 & res$AAB > 8 & res$ABB > 8 & res$AAAB > 8 & res$ABBB > 8) # Double del
		{
			BAlleleFreq <- 0
		}
		else
		{
			BAlleleFreq <- 2
		}
	}
	else if(res2$CNVmeanOrig > 0.1)
	{
		if(res$AAB > 4 & res$AB < 5 & res$ABB > 4) # Dup 3, 10, 3
		{
			BAlleleFreq <- 3
		}
		else if(res$AAAA > 15 & res$AB > 4 & res$AAB > 5 & res$ABB > 4 & res$BBBB > 15) # Double Dup AB = AABB, res$AAAB < 6 & res$ABBB < 6 & res$AB > 5 & res$AAAA > 7 & res$BBBB > 7 & res$AAB > 6 & res$ABB > 6
		{
			BAlleleFreq <- 4
		}
		else
		{
			BAlleleFreq <- 2
		}
	}
	else
	{
		BAlleleFreq <- 2		
	}
	return(BAlleleFreq)
}
