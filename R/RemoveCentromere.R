##' RemoveCentromere: Remove centromeres from the analysis. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RemoveCentromere
##' @param Df: Data frame, Unknown. 
##' @param HG: Human genome version, default = hg19. 
##' @return Unknown.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

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
