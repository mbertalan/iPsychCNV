RemoveCentromere <- function(df, HG="hg18")
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
