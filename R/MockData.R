##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

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
	CNVsSize <- c(25, 50, 100, 150, 300, 450, 600, 900)
	CNVSizeFixed <- sample(CNVsSize, 100, replace=TRUE)
	names(CNVSizeFixed) = 1:100
	
	# Always the same CN
	CNFixed <- sample(c(0,1,1,2,2,3,3,4), 100, replace=TRUE) # CNV Type
	names(CNFixed) = 1:100

	List <- GetMockValues(Type=Type)
	
	Del <- List[["Del"]]
	Dup <- List[["Dup"]] 
	ChrMean <- List[["ChrMean"]]
	ChrMeanProb <- List[["ChrMeanProb"]]
	ChrSD <- List[["ChrSD"]]
	ChrSDProb <- List[["ChrSDProb"]]
	TelomereNoiseSize <- List[["TelomereNoiseSize"]]
	TelomereNoiseEffect <- List[["TelomereNoiseEffect"]]
	BAFs <- List[["BAFs"]]
	BAF_Normal <- List[["BAF_Normal"]]
	BAF_Del <- List[["BAF_Del"]]
	BAF_Dup <- List[["BAF_Dup"]]
	BAF_CN4 <- List[["BAF_CN4"]]
	BAF_CN0 <- List[["BAF_CN0"]]
	BadSNPs <- List[["BadSNPs"]]
	BadSNPIntensity <- List[["BadSNPIntensity"]]
	BadSNPIntensityProb <- List[["BadSNPIntensityProb"]]
	
	
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
			ChrMEAN <- sample(ChrMean[,as.numeric(CHR)], prob=ChrMeanProb[,as.numeric(CHR)], replace=TRUE, size=1)
			X <- sample(ChrMean[,as.numeric(CHR)], prob=ChrMeanProb[,as.numeric(CHR)], replace=TRUE, size=ChrLength)
			
			# Telomere noise
			TelSize1 <- sample(TelomereNoiseSize, 1)
			TelEffect1 <- sample(TelomereNoiseEffect, 1)
			TelSize2 <- sample(TelomereNoiseSize, 1)
			TelEffect2 <- sample(TelomereNoiseEffect, 1)
			# Begining of chromosome
			X[1:TelSize1] <- X[1:TelSize1] + TelEffect1  	
			# End of chrosmoome
			X[(length(X) - TelSize2):length(X)] <- X[(length(X) - TelSize2):length(X)] + TelEffect2
			
			BAF <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=length(X))
			names(BAF) <- SNP.Name
		
			# Adding CNVs		
			NumCNVs <- ((round(length(X)/1000))-1)
			#cat(NumCNVs, "NumCNVs\n")
			DF <- sapply(1:NumCNVs, function(i) # Adding CNVs in the data.
			{
				PositionIndx <- as.numeric(i) * 1000

				# Using fix size for chr position.
				Size <- CNVSizeFixed[i]
				CN <- CNFixed[i]
				
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
				if(CN == 0)
				{
					Impact <- sample(Del, 1)
					BAFCNV <- sample(BAFs, prob=BAF_CN0, replace=TRUE, size=(Size+1))
				}
				if(CN == 4)
				{
					Impact <- sample(Dup, 1)
					BAFCNV <- sample(BAFs, prob=BAF_CN4, replace=TRUE, size=(Size+1))
				}
				if(CN == 2) # Noise data. BAF can not match with LRR
				{
					NoiseLRR <- sample(c(1,3), 1)
					if(NoiseLRR == 1){ Impact <- sample(Del, 1) }
					if(NoiseLRR == 3){ Impact <- sample(Dup, 1) }
					
					if(NoiseLRR == 1)
					{
						NoiseBAF <- sample(c(2,3,4), 1)
						if(NoiseBAF == 3){ BAFCNV <- sample(BAFs, prob=BAF_Dup, replace=TRUE, size=(Size+1)) }
						if(NoiseBAF == 4){ BAFCNV <- sample(BAFs, prob=BAF_CN4, replace=TRUE, size=(Size+1)) }
						if(NoiseBAF == 2){ BAFCNV <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=(Size+1)) }
					}
					else
					{	
						NoiseBAF <- sample(c(2,0,1), 1)
						if(NoiseBAF == 0){ BAFCNV <- sample(BAFs, prob=BAF_CN0, replace=TRUE, size=(Size+1)) }
						if(NoiseBAF == 1){ BAFCNV <- sample(BAFs, prob=BAF_Del, replace=TRUE, size=(Size+1)) }
						if(NoiseBAF == 2){ BAFCNV <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=(Size+1)) }
					}
				}
				
				## Changing GLOBAL VARIABLES ##
				# LRR = X
				X[IndxV] <<- X[IndxV] + Impact
				
				## Changing GLOBAL VARIABLES ##
				BAF[IndxV] <<- BAFCNV
				
				df <- data.frame(Start=Position[PositionIndx], Stop=Position[(PositionIndx+Size)], StartIndx=PositionIndx, StopIndx=(PositionIndx+Size), NumSNPs=Size, Chr=CHR, CNVmean=Impact, CN=CN, sd=SD, ID=FileName, NumCNVs=NumCNVs, ChrMean=ChrMEAN, stringsAsFactors=FALSE)
				return(df)
			})
			df <- MatrixOrList2df(DF)
			#save(df, file="df.RData")
			df2 <- data.frame(SNP.Name=SNP.Name, Chr=rep(CHR, length(X)), Position=Position, Log.R.Ratio=X, B.Allele.Freq=BAF, stringsAsFactors=FALSE)
			return(list(LRR=df2, CNVs=df))
		})
		
		DF <- MatrixOrList2df(tmp["CNVs",])
		LRR <- MatrixOrList2df(tmp["LRR",])
	
		write.table(LRR, sep="\t", quote=FALSE, row.names=FALSE, file=FileName) 
		return(DF)
	})
	
	CNVs <- MatrixOrList2df(All)
	#save(CNVs, file="CNVs.RData")
	CNVs$Length <- CNVs$Stop -  CNVs$Start
	CNVs$CNVID <- 1:nrow(CNVs)
	CNVs$PositionID <- apply(CNVs, 1, function(X){ gsub(" ", "", paste(X["StartIndx"], X["StopIndx"], sep="_", collapse="")) })

	# Adding No-CNVs
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
