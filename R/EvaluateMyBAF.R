##' EvaluateMyBAF: Evaluate B Allele Frequency (BAF). 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title EvaluateMyBAF
##' @param res: Unknown.
##' @param res2. Unknown.
##' @return Classification for BAF.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

EvaluateMyBAF <- function(res, res2)
{
	BAlleleFreq <- 2
	if(res2$CNVmean < 0) # It might be the whole chromosome a CNV
	{
		if(res$AB < 4 & res$ABB < 4 & res$AAB < 4 & res$AAAA > 15  &  & res$BBBB > 15)  # Del
		{
			BAlleleFreq <- 1
		}
		else if(res$AAAA < 25 & res$BBBB < 25 & res$AB > 8 & res$AAB > 8 & res$ABB > 8 & res$AAAB > 8 & res$ABBB > 8) # Double del
		{
			BAlleleFreq <- 0
		}
		else
		{
			BAlleleFreq <- 2
		}
	}
	else if(res2$CNVmean > 0)
	{
		#Centers <- c(0, 0.25, 0.34, 0.5, 0.67, 0.75, 1)
		#names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
		if(res$AAB > 7 & res$AB > 7 & res$ABB > 7 & res$AAAB < 5 & res$ABBB < 5) # Double Dup AB = AABB, res$AAAB < 6 & res$ABBB < 6 & res$AB > 5 & res$AAAA > 7 & res$BBBB > 7 & res$AAB > 6 & res$ABB > 6
		{
			BAlleleFreq <- 4
		}
		else if(res$AAB > 5 & res$ABB > 5 & res$AB < 5 & res$AAAB < 5 & res$ABBB < 5) # Dup 3, 10, 3
		{
			BAlleleFreq <- 3
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
