EvaluateMockResults <- function(MockCNVs, df)
{
	Eval <- apply(MockCNVs, 1, function(X)
	{
		StartM <- as.numeric(X["Start"])
		StopM <- as.numeric(X["Stop"])
		IDM <- X["ID"]
		ChrM <- X["Chr"]
		LengthM <- as.numeric(X["Length"])
		NumSNPsM <- as.numeric(X["NumSNPs"])
		CNVID <- X["CNVID"]

 		res <- subset(df, Chr %in% ChrM & Start < StopM & Stop > StartM & ID %in% IDM)
		if(nrow(res) == 1)
		{
			Found <- TRUE
			StartO <- max(c(res$Start, StartM))
			StopO <- min(c(res$Stop, StopM))
			LengthO <- StopO - StartO

			OverlapLenghM <- ((LengthO/LengthM)*100)
			OverlapSNP <- ((res$NumSNPs/NumSNPsM)*100)
		}
		if(nrow(res) == 0)
		{
			Found <- FALSE
			OverlapLenghM <- 0
			OverlapSNP <-0
		}
		if(nrow(res) > 1)
		{		
			Found <- FALSE
			OverlapLenghM <- 0
			OverlapSNP <-0
		}

		df2 <- data.frame(Found=Found, Overlap.Lengh=OverlapLenghM, Overlap.SNP=OverlapSNP, CNVID, stringsAsFactors=FALSE)
		return(df2)
	})
	Eval <- MatrixOrList2df(Eval)
	df3 <- cbind(Eval, MockCNVs)
	return(df3)
}
