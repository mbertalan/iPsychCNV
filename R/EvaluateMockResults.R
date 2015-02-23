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

		# From df predition
 		res <- subset(df, Chr %in% ChrM & Start < StopM & Stop > StartM & ID %in% IDM)
		CNVID2 <- res$CNVID
		CN2 <- res$CN
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
			CNVID2 <- NA
			CN2 <- NA
		}
		if(nrow(res) > 1)
		{		
			Found <- TRUE
			OverlapLenghM <- 0
			OverlapSNP <-0
			CNVID2 <- NA
			CN2 <- unique(res$CN)[1]
		}

		df2 <- data.frame(Found=Found, Overlap.Length=OverlapLenghM, Overlap.SNP=OverlapSNP, CNVID.Mock=CNVID, CNVID.Pred=CNVID2, CN.Pred=CN2, stringsAsFactors=FALSE)
		return(df2)
	})
	Eval <- MatrixOrList2df(Eval)
	df3 <- cbind(Eval, MockCNVs)
	
	# Generating final df
	df3$Class <- df3$CN
	df3$Class[!df3$Found] <- 2
	
	# CNVs predicted by the method but do not exist.
	FalsePositiveDf <- df[!df$CNVID %in% unique(df3$CNVID.Pred[!is.na(df3$CNVID.Pred)]),]
	if(nrow(FalsePositiveDf) > 0)
	{
		df4 <- apply(FalsePositiveDf, 1, function(Y)
		{
			CN <- -1
			tmpDf <- data.frame(matrix(NA, nrow=1, ncol=ncol(df3)))
			names(tmpDf) <- colnames(df3)
			tmpDf["CNVID.Pred"] <- Y["CNVID"]
			tmpDf["CNVID.Pred"] <- Y["CN"]
			tmpDf["CN"] <- -1 	# Mock does not have CNV there.
			tmpDf["Class"] <- Y["CN"]
			return(tmpDf)
		})
		df4 <- MatrixOrList2df(df4)	
		df3 <- rbind(df3, df4)	
	}
	df3$Class <- as.numeric(df3$Class)
	return(df3)
}
