##' FilterCNVs.V4.: Function to filter predicted Copy Number Variation (CNVs) and avoid a high number of false positive calls. 
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
##'

FilterCNVs.V4 <- function(CNVs = CNVs, MinNumSNPs=20, Sample, ID="Test", verbose=FALSE) #  PathRawData = "~/IBP/CNV/Data/rawData/pilotBroad/"
{	
	CNVID <- rownames(CNVs)
	CNVs$CNVID <- CNVID
	CNV <- Sample
	
	CNVs <- subset(CNVs, NumSNPs > MinNumSNPs)
	
	
	
	AllRes <- apply(CNVs, 1, function(Y) # Loop for CNVs
	{  
		if(verbose){ cat(Y, "\n") }
		CHR <- Y["Chr"]
		CHR <- gsub(" ", "", CHR)
		CNVStart <- as.numeric(Y["Start"]) 
		CNVStop <- as.numeric(Y["Stop"]) 
		NumSNPs <- as.numeric(Y["NumSNPs"])
		Probs <- as.numeric(Y["prob"])
		Size <- CNVStop - CNVStart
		ID <- ID
		if(verbose){ cat(CHR, CNVStart,CNVStop,NumSNPs,Size,ID, Probs, "\n") }

		# Subselection of Data
		tmp <- subset(CNV, Chr %in% CHR) # Whole Chr
		SDChr <- sd(tmp$LRR)
		MeanChr <- mean(tmp$LRR)
		tmp <- tmp[with(tmp, order(Position)), ]
		tmp$PosIndx <- 1:nrow(tmp)
		
		# Getting HMM prob if Probs is null
		if(length(Probs) == 0 || is.na(Probs))
		{
			LRR.mod <- depmix(Log.R.Ratio ~ 1, family = gaussian(), nstates = 3, data = tmp, instart=c(0.1, 0.8, 0.1), respstart=c(-0.45,0,0.3, 0.2,0.2,0.2))
			MyFit2 <- setpars(LRR.mod, getpars(HMM.LRR.fit))
			LRR.probs <- viterbi(MyFit2)
			MeanStates <- apply(LRR.prbs[CNVStart:CNVStop, 2:4], 2, mean) 
			Probs <- sort(MeanStates, decreasing=TRUE)[1] 
			if(verbose){ cat("Probs:", Probs, "\n") }
		}
	
		# GC and LRR
		Chr.SNP.GC <- subset(SNPs.GC, Chr %in% CHR)
		Chr.SNP.GC.LRR <- merge(tmp, Chr.SNP.GC, all=TRUE)
		Chr.SNP.GC.LRR <- subset(Chr.SNP.GC.LRR, !is.na(Log.R.Ratio)) # remove NAs inserted by merged
		Cor.LRR.GC <- cor(Chr.SNP.GC.LRR$GC, Chr.SNP.GC.LRR$LRR, use="p")

		#save(Chr.SNP.GC.LRR, tmp, Chr.SNP.GC, file="Chr.SNP.GC.LRR.RData")
		# Only the CNV region
		tmpRaw <- subset(Chr.SNP.GC.LRR, Position >= CNVStart & Position <= CNVStop)	
		#save(tmpRaw, tmp, file="tmpRaw.RData")
		CNV.LRR.GC <- cor(tmpRaw$GC, tmpRaw$LRR, use="p")
		if(sum(tmpRaw$GC[!is.na(tmpRaw$GC)] < 55) > 20)
		{
			LowGC <- subset(tmpRaw, GC < 55)
			CNVmeanLowGC <- mean(LowGC$LRR)
		}
		else
		{
			CNVmeanLowGC <- mean(tmpRaw$LRR)
		}
			
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
		if(length(Data$High) > 10 & length(Data$CNV) > 10){ CNV2HighPvalue <- t.test(Data$High, Data$CNV)$p.value }else{ CNV2HighPvalue <- 0 }
		if(length(Data$Low) > 10 & length(Data$CNV) > 10) { CNV2LowPvalue <- t.test(Data$Low, Data$CNV)$p.value }else{ CNV2LowPvalue <- 0 }

		ptm.tmp <- proc.time()
		res2 <- GetDataVariables(Data)
		res2$CNVmeanByRef <- res2$CNVmean - MeanChr
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("GetDataVariables time: ", Res.tmp["elapsed"], "\n") }

		# My BAF Classification	
		ptm.tmp <- proc.time()
		res <- ClassNumbers(tmpRaw)
		MyBAF <- EvaluateMyBAF(res, res2)
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("ClassNumbers time: ", Res.tmp["elapsed"], "\n") }
	
		# Defining LogRRatio
		ptm.tmp <- proc.time()
		LogRRatio <- DefiningLogRRatio(res2)
		#if(CNV2HighPvalue < 0.01 || CNV2LowPvalue < 0.01)
		#{
		#	
		#}else{ LogRRatio <- 2 }
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("Define LRR time: ", Res.tmp["elapsed"], "\n") }
		

		# Class by turnpoint: BAlleleFreq by density # Step detection
		ptm.tmp <- proc.time()
		BAFDes <- density(tmpRaw$B.Allele.Freq, adjust = 0.2)
		tp <- turnpoints(BAFDes$y)
		
		# Cleaning Peaks
		tp <- CleaningPeaks(tp)
		SumPeaks <- sum(tp$peaks == TRUE)

		# Defining BAlleleFreq
		dfTmp <- DefineBAFType(SumPeaks)
		BAlleleFreq <- dfTmp$BAlleleFreq
		Class <- dfTmp$Class
		Res.tmp <- proc.time() - ptm.tmp
		if(verbose){ cat("Define BAF time: ", Res.tmp["elapsed"], "\n") }

		# Get genotype Info
		#Genotype <- GetGenotypeInfo(tmpRaw)

		# Add info to res
		#res3 <- cbind(res,res2, Genotype)
		res3 <- cbind(res,res2)
		res4 <- AddInfo2res(res3, CNV2HighPvalue, CNV2LowPvalue, Class, BAlleleFreq, MyBAF, LogRRatio, SumPeaks, SDChr, MeanChr)
		res4$Cor.LRR.GC <- Cor.LRR.GC  
		res4$CNVmeanLowGC <- CNVmeanLowGC
		if(is.null(Y["prob"])){ res4$prob <- Probs }
		return(res4)
	})
	df <- do.call(rbind, AllRes)
	tmp2 <- cbind(CNVs, df) # Combining position variables with filter ones.
	#save(tmp2, file="tmp2.RData")
	# Define if the CNV is Good or Bad
	#Type <- DefineCNVType(tmp2)
	#tmp2$Type <- Type
	Class <- DefineCNVClass(tmp2)
	tmp2$Class <- Class
	df <- tmp2
	df$Source <- rep("iPsychCNV", nrow(df))
	df$CN <- df$Class
	df$CN[df$CN %in% "Del"] <- "1"
	df$CN[df$CN %in% "Normal"] <- "2"
	df$CN[df$CN %in% "Dup"] <- "3"
	df$CN[df$CN %in% "DoubleDel"] <- "0"
	df$CN[df$CN %in% "DoubleDup"] <- "4"
	df$CN <- as.numeric(df$CN)
	df <- df[, !colnames(df) %in% "Class"]
	df$ID <- ID
	return(df)
}

