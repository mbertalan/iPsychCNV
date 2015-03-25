##' Filter CNVs: function to filter predicted CNVs and avoid high number of false positives. 
##'
##' The function receives a data frame with CNV information, ex: Chr, Start position, Stop position and Sample ID. Returns a data frame with CNVs and a Type difining if is bad or good.
##' @title iPsychCNVs
##' @return Data frame with CNVs and classification.
##' @author Marcelo Bertalan
##' @export

FilterCNVs.V4 <- function(CNVs = CNVs, MinNumSNPs, CNV, ID) #  PathRawData = "~/IBP/CNV/Data/rawData/pilotBroad/"
{	
	CNVID <- rownames(CNVs)
	CNVs$CNVID <- CNVID
	
	AllRes <- apply(CNVs, 1, function(Y) # Loop for CNVs
	{  
		CHR <- as.numeric(Y["Chr"])
		CNVStart <- as.numeric(Y["Start"]) 
		CNVStop <- as.numeric(Y["Stop"]) 
		NumSNPs <- as.numeric(Y["NumSNPs"])
		Size <- CNVStop - CNVStart
		ID <- ID

		# Subselection of Data
		tmp <- subset(CNV, Chr %in% CHR) # Whole Chr
		SDChr <- sd(tmp$LRR)
		MeanChr <- mean(tmp$LRR)
		tmp <- tmp[with(tmp, order(Position)), ] 
		tmp$PosIndx <- 1:nrow(tmp)

		tmpRaw <- subset(tmp, Position >= CNVStart & Position <= CNVStop)	# Only the CNV region
			
		
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
		if(length(Data$High) > 20 & length(Data$CNV) > 20){ CNV2HighPvalue <- t.test(Data$High, Data$CNV)$p.value }else{ CNV2HighPvalue <- 0 }
		if(length(Data$Low) > 20 & length(Data$CNV) > 20) { CNV2LowPvalue <- t.test(Data$Low, Data$CNV)$p.value }else{ CNV2LowPvalue <- 0 }

		#cat("Get Variables", " ", proc.time(), "\n")
		res2 <- GetDataVariables(Data)

		# My BAF Classification	
		res <- ClassNumbers(tmpRaw$B.Allele.Freq)
		MyBAF <- EvaluateMyBAF(res, res2$CNVmean)
	
		# Defining LogRRatio
		if(CNV2HighPvalue < 0.01 || CNV2LowPvalue < 0.01)
		{
			LogRRatio <- DefiningLogRRatio(res2)
		}else{ LogRRatio <- "Undefined" }
		

		# Class by turnpoint: BAlleleFreq by density # Step detection
		#cat("Turn points", " ", proc.time(), " \n")
		BAFDes <- density(tmpRaw$B.Allele.Freq,adjust=0.5)
		tp <- turnpoints(BAFDes$y)
		
		# Cleaning Peaks
		tp <- CleaningPeaks(tp)
		SumPeaks <- sum(tp$peaks == TRUE)

		# Defining BAlleleFreq
		dfTmp <- DefineBAFType(SumPeaks)
		BAlleleFreq <- dfTmp$BAlleleFreq
		Class <- dfTmp$Class

		# Add info to res
		res3 <- cbind(res,res2)
		res4 <- AddInfo2res(res3, CNV2HighPvalue, CNV2LowPvalue, Class, BAlleleFreq, MyBAF, LogRRatio, SumPeaks, SDChr, MeanChr)
		return(res4)
	})
	df <- do.call(rbind, AllRes)
	tmp2 <- cbind(CNVs, df) # Combining position variables with filter ones.

	
	# Define if the CNV is Good or Bad
	Type <- DefineCNVType(tmp2)
	Class <- DefineCNVClass(tmp2)
	tmp2$Type <- Type
	tmp2$Class <- Class
	return(tmp2)
}

