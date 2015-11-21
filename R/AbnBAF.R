##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##' Specifically designed to reduce false positive CNVs and handle data from amplified DNA on dried blood spots.
##'
##' AbnBAF

##' @title AbnBAF
##' @param res
##' @param LRR
##' @param chrY.Perce

##' @return BAlleleFreq: Estimate BAF copy number. 
##' @author Marcelo Bertalan
##' @export
##' @examples
##' mockCNV <- MockData(N=5, Type="Blood", Cores=1)
##' cnvs <- iPsychCNV(PathRawData=".", Cores=1, Pattern="^MockSample*", Skip=0)


AbnBAF <- function(res, LRR, chrY.Perc)
{
	BAlleleFreq <- -1

	if(res$AAAA >= 15  & res$AB <= 5 & res$BBBB >= 15 & res$ABB <= 5 & res$AAB <= 5 & LRR > -1 & LRR < 0.2 & chrY.Perc <= 15)  # Del
	{
			BAlleleFreq <- 1
	}
	else if(res$AAAA <= 25 & res$BBBB <= 25 & res$AB >= 8 & res$AAB >= 5 & res$ABB >= 5 & res$AAAB <= 25 & res$ABBB <= 25 & LRR <= -0.6 & chrY.Perc > 30) # Double del
	{
		BAlleleFreq <- 0
	}
	else if(res$AAAA >= 15 & res$BBBB >= 15 & res$AB >= 5 & LRR > -0.2 & chrY.Perc >= 9)
	{
		BAlleleFreq <- 2
	}
	else if(res$AAB >= 4 & res$AB <= 5 & res$ABB >= 4 & LRR > 0.1  & chrY.Perc >= 5) 
	{
		BAlleleFreq <- 3
	}
	else if(res$AAAA >= 15 & res$AB >= 4 & res$AAB >= 5 & res$ABB >= 4 & res$BBBB >= 15 & LRR > 0.2 & chrY.Perc >= 5) 
	{
		BAlleleFreq <- 4
	}
	else 
	{
		# define by LRR only
		if(LRR > -1 & LRR <= 0.1 & chrY.Perc <= 15)
		{
			BAlleleFreq <- 1
		}
		else if(LRR >= 0.1  & chrY.Perc >= 15)
		{
			BAlleleFreq <- 2
		}
		else if(LRR >= 0.1  & chrY.Perc >= 15)
		{
			BAlleleFreq <- 1
		}
		else if(LRR <= -0.6 & chrY.Perc > 15)
		{
			BAlleleFreq <- 0
		}
	}
	
	if(is.na(BAlleleFreq)) # If still NA Likely to be zero or 1 from Y or mt
	{
		if(res$AAAA >= 15  & res$AB <= 5 & res$BBBB >= 15 & res$ABB <= 5 & res$AAB <= 5 & chrY.Perc <= 15) # Double del
		{
			BAlleleFreq <- 1
		}
	}
	
	return(BAlleleFreq)
}
