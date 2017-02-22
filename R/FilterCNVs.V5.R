##' FilterCNVs.V5.: Function to filter predicted Copy Number Variation (CNVs) and avoid a high number of false positive calls. 
##'
##' The function receives a data frame with CNV information, ex: chr., start position, stop position, and sample ID. 
##' @title FilterCNVs.V4.
##' @param CNVs: Data frame with CNVs. Unknown?
##' @param MinNumSNPs: Minimum number of SNPs per CNV, default = 20.
##' @param Sample: Unknown.
##' @param ID: Unknown.
##' @param Verbose: Unknown, default = FALSE.
##' @return Data frame with CNVs and classification. 
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.

FilterCNVs.V5 <- function(CNVs = CNVs, MinNumSNPs=20, Sample, ID="Test", verbose=FALSE) #  PathRawData = "~/IBP/CNV/Data/rawData/pilotBroad/"
{	
	suppressPackageStartupMessages(library("depmixS4"))
	library(fpc)
	
	CNVID <- rownames(CNVs)
	CNVs$CNVID <- CNVID
	CNV <- Sample
	CNVs <- subset(CNVs, NumSNPs > MinNumSNPs)
	ptm.tmp <- proc.time()
	
	AllRes <- apply(CNVs, 1, function(Y) # Loop for CNVs
	{  
		if(verbose){ cat(Y, "\n") }
		CHR <- Y["Chr"]
		CHR <- gsub(" ", "", CHR)
		CNVStart <- as.numeric(Y["Start"]) 
		CNVStop <- as.numeric(Y["Stop"])
		IndxStart <- as.numeric(Y["StartIndx"])
		IndxStop <- as.numeric(Y["StopIndx"])
		NumSNPs <- as.numeric(Y["NumSNPs"])
		Size <- CNVStop - CNVStart
		ID <- ID
		if(verbose){ cat("# Start #", CHR, CNVStart,CNVStop,NumSNPs,Size,ID, "\n") }

		# Subselection of Data
		tmp <- subset(CNV, Chr %in% CHR) # Whole Chr
		SDChr <- sd(tmp$LRR)
		MeanChr <- mean(tmp$LRR)
		tmp <- tmp[with(tmp, order(Position)), ]
		tmp$PosIndx <- 1:nrow(tmp)
		tmpRaw <- subset(tmp, Position >= CNVStart & Position <= CNVStop)	
		#save(tmpRaw, tmp, file="tmpRaw.RData")
			

		# Before and after CNV
		IndxStart <- tmpRaw$PosIndx[1] - NumSNPs
		if(IndxStart < 0){ IndxStart <- 1; LowMean <- 0 }
		IndxStop <- tmpRaw$PosIndx[length(tmpRaw$PosIndx)] + NumSNPs
		if(IndxStop > nrow(tmp)){  IndxStop <- nrow(tmp); HighMean <- 0 }

		CNVStartIndx <- tmpRaw$PosIndx[1]
		CNVStopIndx <- tmpRaw$PosIndx[length(tmpRaw$PosIndx)]

		Low <- subset(tmp, PosIndx >= IndxStart &  PosIndx <= CNVStartIndx)	# Selecting data before CNV
		High <- subset(tmp, PosIndx >= CNVStopIndx &  PosIndx <= IndxStop)	# Selecting data after CNV
		
		# Creating a list with variables
		Data <- list(Low=Low$Log.R.Ratio, CNV=tmpRaw$Log.R.Ratio, High=High$Log.R.Ratio, LRR=tmpRaw$LRR) # LRR is original corrected for mean = 0.
		#save(Data, file="Data.RData")
		if(length(Data$High) > 10 & length(Data$CNV) > 10){ CNV2HighPvalue <- t.test(Data$High, Data$CNV)$p.value }else{ CNV2HighPvalue <- 0 }
		if(length(Data$Low) > 10 & length(Data$CNV) > 10) { CNV2LowPvalue <- t.test(Data$Low, Data$CNV)$p.value }else{ CNV2LowPvalue <- 0 }

		# Get Info
		ptm.tmp <- proc.time()
		res2 <- GetDataVariables(Data)
		res2$CNVmeanByRef <- res2$CNVmean - MeanChr
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("GetDataVariables time: ", Res.tmp["elapsed"], "\n") }

		# BAF classification by Partitioning Around Medoids
		ptm.tmp <- proc.time()
			# removing BAF == 0 and == 1 to increase speed if CNV has more than 200 SNPs
		if(nrow(tmpRaw) > 200){ BAF <- tmpRaw$B.Allele.Freq[tmpRaw$B.Allele.Freq != 0 & tmpRaw$B.Allele.Freq != 1] }
		else{ BAF <- tmpRaw$B.Allele.Freq }

		M <- sapply(BAF, function(I){ I - BAF })
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("sapply time: ", Res.tmp["elapsed"], "\n") }
		
		ptm.tmp <- proc.time()
		pamk.best <- pamk(M)
		PamBAF <- pamk.best$nc
		BAlleleFreq <- as.numeric(DefineBAFType(PamBAF))
		if(is.na(BAlleleFreq)){ BAlleleFreq <- 2 }
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("PAM time: ", Res.tmp["elapsed"], "\n") }
		
		# My BAF Classification	
		ptm.tmp <- proc.time()
		res <- ClassNumbers(tmpRaw)
		MyBAF <- EvaluateMyBAF(res, res2)
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("ClassNumbers time: ", Res.tmp["elapsed"], "\n") }
	
		# Defining LogRRatio
		ptm.tmp <- proc.time()
		LogRRatio <- DefiningLogRRatio(res2)
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("Define LRR time: ", Res.tmp["elapsed"], "\n") }
		
		res3 <- cbind(res,res2)
		res3$BAlleleFreq <- BAlleleFreq
		res3$MyBAF <- MyBAF
		res3$LogRRatio <- LogRRatio 
		res3$SDChr <- SDChr
		res3$MeanChr <- MeanChr
		res3$PamBAF <- BAlleleFreq

		return(res3)
	})
	Res.tmp <- proc.time() - ptm.tmp
	if(verbose){ cat("Done: ", Res.tmp["elapsed"], "\n") }
	system.time(df <- MatrixOrList2df(AllRes))
	tmp2 <- cbind(CNVs, df) # Combining position variables with filter ones.
	
	Class <- DefineCNVClass(tmp2)
	tmp2$CN <- as.numeric(Class)
	tmp2$Source <- "iPsychCNV"
	tmp2$ID <- ID
	
	return(tmp2)
}
