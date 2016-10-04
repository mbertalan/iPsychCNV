AddNoCNVRegions <- function(CNVs, CNV)
{
	CNVs$PositionID <- apply(CNVs, 1, function(X){ gsub(" ", "", paste(X["StartIndx"], X["StopIndx"], sep="_", collapse="")) })

	tmp3 <- sapply(unique(CNVs$ID), function(X)
	{
		tmp <- sapply(unique(CNVs$Chr), function(Y)
		{
			# selecting the CNVs
			subCNVs <- subset(CNVs, ID %in% X & Chr %in% Y)
			# Selecting data info from chip
			subCNV <- subset(CNV, Chr %in% Y)
			
			indx <- sort(c(1, subCNVs$StartIndx, subCNVs$StopIndx, nrow(subCNV)))
			Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
			Info <- t(Info)
			Df <- as.data.frame(Info)
			names(Df) <- c("StartIndx", "StopIndx")
			PositionID <- apply(Df, 1, function(X){ gsub(" ", "", paste(X[1], X[2], sep="_", collapse="")) })
			Df$PositionID <- PositionID
			Df2 <- Df[!Df$PositionID %in% subCNVs$PositionID,] # Only the non-CNVs
			subCNV <- subset(CNV, Chr %in% Y)
			subCNV <- subCNV[order(subCNV$Position),]
			Df2$Start <- subCNV$Position[Df2$StartIndx]
			Df2$Stop <- subCNV$Position[Df2$StopIndx]
			Df2$NumSNPs <- Df2$StopIndx - Df2$StartIndx
			Df2$Chr <- rep(Y, nrow(Df2))
			Df2$CN <- rep(2, nrow(Df2))
			Df2$ID <- rep(X, nrow(Df2))
			Df2$Length <- Df2$Stop - Df2$Start
			#Df2 <- Df2[, colnames(subCNVs)]
			return(Df2)
		})
		tmp2 <- MatrixOrList2df(tmp)
		tmp2 <- tmp2[, !colnames(tmp2) %in% ".id"]
		return(tmp2)
	})
	tmp4 <- MatrixOrList2df(tmp3)
	tmp5 <- rbind(tmp4, CNVs[, colnames(tmp4)])
	tmp5$CNVID <- 1:nrow(tmp5)
	return(tmp5)
}
