EvaluateMockResults <- function(MockCNVs, df)
{
	library(pROC)
	Eval <- apply(MockCNVs, 1, function(X)
	{
		StartM <- as.numeric(X["Start"])
		StopM <- as.numeric(X["Stop"])
		IDM <- X["ID"]
		ChrM <- X["Chr"]
		LengthM <- as.numeric(X["Length"])
		NumSNPsM <- as.numeric(X["NumSNPs"])
		CNVID <- gsub(" ", "", X["CNVID"])
		CNM <- as.numeric(X["CN"])

		if(CNM == 2){ CNV.Present=0 }else{ CNV.Present=1 }
		# From df predition
 		res <- subset(df, Chr %in% ChrM & Start < StopM & Stop > StartM & ID %in% IDM)
		CNVID2 <- res$CNVID
		CN2 <- res$CN
		if(nrow(res) == 1)
		{
			if(CN2 == CNM){ Found <- 1 }else{ Found <- 0}
			StartO <- max(c(res$Start, StartM))
			StopO <- min(c(res$Stop, StopM))
			LengthO <- StopO - StartO

			OverlapLenghM <- ((LengthO/LengthM)*100)
			OverlapSNP <- ((res$NumSNPs/NumSNPsM)*100)
		}
		if(nrow(res) == 0)
		{
			Found <- 0
			OverlapLenghM <- 0
			OverlapSNP <-0
			CNVID2 <- NA
			CN2 <- 2
		}
		if(nrow(res) > 1)
		{
			CN2 <- CN2[1]
			if(CN2 == CNM){ Found <- 1 }else{ Found <- 0}	
			OverlapLenghM <- 0
			OverlapSNP <-0
			CNVID2 <- CNVID[1]
		}

		df2 <- data.frame(CNV.Present=CNV.Present, CNV.Found=Found, Overlap.Length=OverlapLenghM, Overlap.SNP=OverlapSNP, CNVID.Mock=CNVID, CNVID.Pred=CNVID2, CN.Pred=CN2, stringsAsFactors=FALSE)
		return(df2)
	})
	Eval <- MatrixOrList2df(Eval)
	df3 <- cbind(Eval, MockCNVs)
	return(df3)
}
