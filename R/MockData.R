MockData <- function(N=1, Wave1=FALSE, Type="Blood", Cores=1) # Type: Blood or PKU (Perfect or noise)
{
	# Use Wave1 PFB ? Wave1PFB comes with the package.
	if(!Wave1)
	{
		WaveTmp <- rep(0.5, length(Wave1PFB))
		names(WaveTmp) <- names(Wave1PFB)
		Wave1PFB <- WaveTmp
	}
		
	# CNV Info: Using always the same position. Rdata from the package.
	CNVsSize <- c(10, 30, 50, 100, 150, 200, 300, 400)
	#CNVSizeFixed <- sample(CNVsSize, 50, replace=TRUE)
	#names(CNVSizeFixed) = 1:50

	List <- GetMockValues(Type=Type)
	
	Del <- List[["Del"]]
	Dup <- List[["Dup"]] 
	ChrMean <- List[["ChrMean"]]
	ChrSD <- List[["ChrSD"]]
	ChrSDProb <- List[["ChrSDProb"]]
	TelomereNoiseSize <- List[["TelomereNoiseSize"]]
	TelomereNoiseEffect <- List[["TelomereNoiseEffect"]]
	BAFs <- List[["BAFs"]]
	BAF_Normal <- List[["BAF_Normal"]]
	BAF_Del <- List[["BAF_Del"]]
	BAF_Dup <- List[["BAF_Dup"]]
	BadSNPs <- List[["BadSNPs"]]
	BadSNPIntensity <- List[["BadSNPIntensity"]]
	BadSNPIntensityProb <- List[["BadSNPIntensityProb"]]
	ChrMeanProb <- List[["ChrMeanProb"]]
	
	if(Type %in% "Blood")
	{
		BAF_Prob_By_Chr <- ((55*GC_MeanByChr)/40)/100	
	}
	else
	{
		BAF_Prob_By_Chr <- (((65*GC_MeanByChr)/40)+(GC_MeanByChr-45))/100
	}
	GC_ByChr <- GC_MeanByChr + ((GC_MeanByChr-median(GC_MeanByChr))*4)
	GC_Effect <- (GC_ByChr/median(GC_ByChr))^2
	
	suppressPackageStartupMessages(library(parallel))

	All <- mclapply(1:N, mc.cores=Cores, mc.preschedule = FALSE, function(SampleNum) 
	{
		FileName <- paste("MockSample_", SampleNum, ".tab", sep="", collapse="")
	
		tmp <- sapply(unique(CNV$Chr), function(CHR) # Chromosome loop
		{
			# CNV is a RData from the package with the Psycho chip information. 
			# It has SNP position and names. It is used to guide the mock data.
			subCNV <- subset(CNV, Chr %in% CHR) 
			subCNV <- subCNV[order(subCNV$Position),]
			Position <- subCNV$Position
			SNP.Name <- subCNV$SNP.Name
	
			ChrLength <- nrow(subCNV)
			SD=sample(ChrSD, 1, prob=ChrSDProb) # chr sd
			ChrMEAN <- sample(ChrMean, prob=ChrMeanProb, replace=TRUE, size=1)
			#X <- rnorm(ChrLength, sd=SD, mean=ChrMEAN)
			X <- sample(ChrMean, prob=ChrMeanProb, replace=TRUE, size=ChrLength)
			X <- as.numeric(X)
			
			# Change BAF to simulate chromosome differences
			Tmp_BAF_Prob <- BAF_Normal
			
			#Change BAF frequency (High on 0 or High on 1). It seems High on 1 give dupliations on penncnv.
			BAF_Change <- sample(c(1,2), 1)
			if(BAF_Change == 1)
			{ 
				BAF_Prob_Value_BBBB <- BAF_Prob_By_Chr[CHR]
				BAF_Prob_Value_AAAA <- 1 - BAF_Prob_Value_BBBB
			}
			else
			{
				BAF_Prob_Value_AAAA <- BAF_Prob_By_Chr[CHR]
				BAF_Prob_Value_BBBB <- 1 - BAF_Prob_Value_AAAA
			}
			Tmp_BAF_Prob[100:101] <- BAF_Prob_Value_BBBB
			Tmp_BAF_Prob[98:99] <- (BAF_Prob_Value_BBBB * 3/4)
			Tmp_BAF_Prob[1:2] <- BAF_Prob_Value_AAAA
			Tmp_BAF_Prob[3:4] <- (BAF_Prob_Value_AAAA * 3/4)
			
			BAF <- sample(BAFs, prob=Tmp_BAF_Prob, replace=TRUE, size=length(X))
			names(BAF) <- SNP.Name
			
			# Adding Psych Chip PFB to low freq SNPs. Using fix position
			# Wave1PFB: data from package. Pop frequency estimated by Wave1. 
			IndxBAF1 <- names(BAF) %in% names(Wave1PFB)[Wave1PFB > 0.9]
			if(sum(IndxBAF1) > 1)
			{
				BAF[IndxBAF1] <- rnorm(sum(IndxBAF1), mean=0.97, sd=0.01)
			}
			
			IndxBAF0 <- names(BAF) %in% names(Wave1PFB)[Wave1PFB < 0.1]
			if(sum(IndxBAF0) > 1)
			{
				BAF[IndxBAF0] <- rnorm(sum(IndxBAF0), mean=0.01, sd=0.01)
			}
			
			# Adding random noise
			t  <- 1:length(X)
			ssp <- spectrum(X, plot=FALSE)  
			per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
			reslm <- lm(X ~ sin(2*pi/per*t)+cos(2*pi/per*t))
			if(Type %in% "PKU")
			{
				X <- X - abs((fitted(reslm)*20))
			}
			else
			{
				X <- X - abs((fitted(reslm)*2))
			}

			# Add Telomere noise
			NTelomereSize <- sample(TelomereNoiseSize, 1)
			TeloEffect <- sample(TelomereNoiseEffect, 1) 
			if(Type %in% "PKU"){ TeloEffect <- TeloEffect * GC_Effect[CHR] }else{ TeloEffect <- (TeloEffect/2) * GC_Effect[CHR] }  # If PKU the telomereNoise is double.
			X[1:NTelomereSize] <- X[1:NTelomereSize] + TeloEffect
			X[(length(X) - NTelomereSize):length(X)] <- X[(length(X) - NTelomereSize):length(X)] + TeloEffect
			

			# Adding bad SNPs (generally because of GC and LCR)
			TotalNumberofBadSNPs <- round(length(X)*BadSNPs[CHR])
			TotalNumberofBadSNP <- round(TotalNumberofBadSNP/100) + 1
			BadSNPsIndx <- sample(1:length(X), TotalNumberofBadSNPs)
			BadSNPsIndx <- as.vector(sapply(BadSNPsIndx, function(I){ c((I-40):(I+40)) }))
			BadSNPsIndx[BadSNPsIndx < 0] <- 1
			BadSNPsIndx[BadSNPsIndx > length(X)] <- length(X)
			NoiseSNP <- sample(BadSNPIntensity, prob=BadSNPIntensityProb, 1)
			X[BadSNPsIndx] <- X[BadSNPsIndx] - abs(rnorm(length(BadSNPsIndx), sd=(SD*1.5), mean=NoiseSNP))

			
			
		
			# Adding CNVs		
			NumCNVs <- ((round(length(X)/2000))-2)
			
			DF <- sapply(1:NumCNVs, function(i) # Adding CNVs in the data.
			{
				CN <- sample(c(1,3), 1) # CNV Type
				PositionIndx <- as.numeric(i) * 2000
				#Size <- sample(CNVsSize, 1) # CNV Size
				# Using fix size for chr position.
				Size <- CNVSizeFixed[i]
				
				# CNV position
				IndxV <- PositionIndx:(PositionIndx+Size)
				
				if(CN == 1)
				{
					Impact <- sample(Del, 1)
					BAFCNV <- sample(BAFs, prob=BAF_Del, replace=TRUE, size=(Size+1))
				}
				if(CN == 3)
				{
					Impact <- sample(Dup, 1)
					BAFCNV <- sample(BAFs, prob=BAF_Dup, replace=TRUE, size=(Size+1))
				}
				
				## Changing GLOBAL VARIABLES ##
				# LRR = X
				X[IndxV] <<- X[IndxV] + Impact
	
				# BAF, Change BAF but keep SNPs with low heterozygosity.
				names(BAFCNV) <- IndxV
				NoChangeIndx <- c(which(BAF[IndxV] > 0.9), which(BAF[IndxV] < 0.1))
				NewIndx <- IndxV[(NoChangeIndx*-1)]
				BAF[NewIndx] <<- BAFCNV[as.character(NewIndx)]
				
				## Changing GLOBAL VARIABLES ##
				BAF[BAF > 1] <<- 1
				BAF[BAF < 0] <<- 0
				
				df <- data.frame(Start=Position[PositionIndx], Stop=Position[(PositionIndx+Size)], StartIndx=PositionIndx, StopIndx=(PositionIndx+Size), NumSNPs=Size, Chr=CHR, CNVmean=Impact, CN=CN, sd=SD, ID=FileName, NoiseSNP=NoiseSNP, BadSNPs=TotalNumberofBadSNPs, NumCNVs=NumCNVs, ChrMean=ChrMEAN, stringsAsFactors=FALSE)
				return(df)
			})
			df <- MatrixOrList2df(DF)
			df2 <- data.frame(SNP.Name=SNP.Name, Chr=rep(CHR, length(X)), Position=Position, Log.R.Ratio=X, B.Allele.Freq=BAF, stringsAsFactors=FALSE)
			return(list(LRR=df2, CNVs=df))
		})
		DF <- MatrixOrList2df(tmp["CNVs",])
		LRR <- MatrixOrList2df(tmp["LRR",])
	
		write.table(LRR, sep="\t", quote=FALSE, row.names=FALSE, file=FileName) 
		return(DF)
	})
	CNVs <- MatrixOrList2df(All)
	CNVs$Length <- CNVs$Stop -  CNVs$Start
	CNVs$CNVID <- 1:nrow(CNVs)
	CNVs$PositionID <- apply(CNVs, 1, function(X){ gsub(" ", "", paste(X["StartIndx"], X["StopIndx"], sep="_", collapse="")) })

	# Adding No-CNVs
	tmp3 <- sapply(unique(CNVs$ID), function(X)
	{
		tmp <- sapply(unique(CNVs$Chr), function(Y)
		{
			subCNVs <- subset(CNVs, ID %in% X & Chr %in% Y)
			indx <- sort(c(1, subCNVs$StartIndx, subCNVs$StopIndx))
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
			Df2$CNVmean <- rep(0, nrow(Df2))
			Df2$CN <- rep(2, nrow(Df2))
			Df2$sd <- rep(0.2, nrow(Df2))
			Df2$ID <- rep(X, nrow(Df2))
			Df2$NoiseSNP <- rep(unique(subCNVs$NoiseSNP), nrow(Df2))
			Df2$BadSNPs <- rep(unique(subCNVs$BadSNPs), nrow(Df2))
			Df2$NumCNVs <- rep(unique(subCNVs$NumCNVs), nrow(Df2))
			Df2$ChrMean <- rep(unique(subCNVs$ChrMean), nrow(Df2))
			Df2$Length <- Df2$Stop - Df2$Start
			Df2$CNVID <- 1:nrow(Df2)    
			Df2 <- Df2[, colnames(subCNVs)]
			return(Df2)
		})
		tmp2 <- MatrixOrList2df(tmp)
		tmp2 <- tmp2[, !colnames(tmp2) %in% ".id"]
		return(tmp2)
	})
	tmp3 <- MatrixOrList2df(tmp3)
	tmp3 <- tmp3[, !colnames(tmp3) %in% ".id"]

	tmp4 <- rbind(tmp3, CNVs)
	tmp4 <- tmp4[order(tmp4$ID, tmp4$Chr, tmp4$Start),]
	return(tmp4)
}	
