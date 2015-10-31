GetIndxPositionFromChips <- function(CNVs, Sample)
{
	Res <- apply(CNVs, 1, function(X)
	{
		#cat(X, "\n")
		StartM <- as.numeric(X["Start"]) # CNV Start
		StopM <- as.numeric(X["Stop"])
		CHR <- X["Chr"]
		#cat("Sample rows:", nrow(Sample), "\n")
		CHR <- gsub(" ", "", CHR)	
		subSample <- subset(Sample, Chr %in% CHR) 			# Subset by Chr
		subSample <- subSample[with(subSample, order(Position)),]	
		#cat("SubSample rows:", nrow(subSample), "\n")
		# Order
		StartIndx <- which.min(abs(subSample$Position - StartM))
		StopIndx <- which.min(abs(subSample$Position - StopM))
		NumSNPs <- StopIndx - StartIndx
		#cat(StartIndx, StopIndx, "\n")
		df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx, NumSNPs=NumSNPs)
		return(df)
	})
	#save(Res, file="Res.RData")
	df <- MatrixOrList2df(Res)
	df2 <- cbind(CNVs, df)	
	return(df2)
}
