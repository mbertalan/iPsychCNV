GetIndxPositionFromChips <- function(CNVs, Sample)
{
	Res2 <- sapply(unique(CNVs$Chr), function(y)
	{
		CHR <- gsub(" ", "", y)
		df2 <- subset(CNVs, Chr %in% CHR)	
		subSample <- subset(Sample, Chr %in% CHR) 			# Subset by Chr
		subSample <- subSample[with(subSample, order(Position)),]	
	
		Res <- apply(df2, 1, function(X)
		{
			StartM <- as.numeric(X["Start"]) # CNV Start
			StopM <- as.numeric(X["Stop"])
			# Order
			StartIndx <- which.min(abs(subSample$Position - StartM))
			StopIndx <- which.min(abs(subSample$Position - StopM))
			df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx)
			return(df)
		})
		df <- MatrixOrList2df(tmp=Res)
		return(df)
	})
	df <- MatrixOrList2df(tmp=Res2)
	CNVs <- cbind(CNVs, df)
	CNVs$NumSNPs <- CNVs$StopIndx - CNVs$StartIndx
	return(CNVs)
}
