GetIndxPositionFromChips <- function(CNVs, Sample)
{
	Res2 <- sapply(unique(CNVs$Chr), function(y)
	{
		CHR <- y
		CHR <- gsub(" ", "", CHR)
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
			NumSNPs <- StopIndx - StartIndx
			df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx, NumSNPs=NumSNPs)
			df2 <- cbind(X, df)
			return(df2)
		})
		df <- MatrixOrList2df(tmp=Res)
		return(df)
	})
	df <- MatrixOrList2df(tmp=Res2)
	return(df)
}
