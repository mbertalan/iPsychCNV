##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

RemoveCentromere <- function(df, HG="hg19")
{
	Cent <- subset(Cent, hg %in% HG)
	tmp <- sapply(unique(df$Chr), function(X)
	{
		Cent.tmp <- subset(Cent, Chr %in% X)
		StartC <- Cent.tmp$Start
		StopC <- Cent.tmp$Stop
		
		# selecting the chr
		df.tmp <- subset(df, Chr %in% X)
		
		# removing CNVs that overlap centromere
		df.tmp <- subset(df.tmp, Stop < StartC | Start > StopC)
		
		return(df.tmp)
	})
	tmp2 <- MatrixOrList2df(tmp)
	return(tmp2)
}
