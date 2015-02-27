MockData <- function(N=1)
{
	# CNV Info
	CNVsSize <- c(30, 50, 100, 150, 200, 300, 400, 500)
	CNVSizeFixed <- sample(CNVsSize, 50, replace=TRUE)
	names(CNVSizeFixed) = 1:50
	Del <- seq(from=-0.15, to=-0.45, by=-0.05)
	Dup <- seq(from=0.15, to=0.45, by=0.05)
	CNVMean <- seq(from=0.2, to=-0.4, by=-0.1)
	
	# Telomere noise
	TelomereNoiseSize <- seq(from=100, to=500, by=100)
	TelomereNoiseEffect <- seq(from=-0.1, to=-1, by=-0.1)
	
	# BAF
	BAFs <- seq(from=0, to=1, by=0.01) # 21
	BAF_Basic <- rep(0.02, 101)
	names(BAF_Basic) <- BAFs
	BAF_Basic[10:90] <- seq(from=0.18, to=0.22, by=0.0005)
	
	# BAF normal prob	
	BAF_Normal <- BAF_Basic
	BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + 0.2
	BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + 0.1
	BAF_Normal[c(8:9)] <- BAF_Normal[c(8:9)] - 0.02
	
	BAF_Normal[c(80:85)] <- BAF_Normal[c(80:85)] + 0.01
	BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + 0.20
	BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + 0.5	
	
	# BAF Del prob
	BAF_Del <- BAF_Basic
	BAF_Del[c(1:2)] <- BAF_Normal[c(1:2)] + 0.23
	BAF_Del[c(98:99)] <- BAF_Del[c(98:99)] + 0.2
	BAF_Del[c(100:101)] <- BAF_Del[c(100:101)] + 0.45
	
	# BAF Dup prob
	BAF_Dup <-  BAF_Basic
	BAF_Dup[25:35] <- BAF_Dup[25:35] + 0.18 # 0.25 0.30 0.35, 
	BAF_Dup[70:80] <- BAF_Dup[70:80] + 0.18 # 0.7 0.75 0.8
	
	# BadSNPs
	BadSNPs <- c(0.01, 0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.09,0.08,0.01,0.15,0.05,0.05,0.10)
	names(BadSNPs) <- 1:22
	BadSNPIntensity <- seq(from=0.2, to=-4, by=-0.1)
	BadSNPIntensityProb <- seq(from=0.53, to=0.11, by=-0.01)
	
	All <- sapply(1:N, function(SampleNum) # File loop
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
			SD=sample(seq(from=0.1, to=0.5, by=0.1), 1, prob=c(0.2,0.3,0.3,0.2,0.1)) # chr sd
			MEAN <- sample(CNVMean, 1)
			X <- rnorm(ChrLength, sd=SD, mean=MEAN)
			BAF <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=length(X))
			names(BAF) <- SNP.Name
			
			# Adding Psych Chip PFB to low freq SNPs. Using fix position
			IndxBAF1 <- names(BAF) %in% names(Wave1PFB)[Wave1PFB > 0.95]
			if(sum(IndxBAF1) > 1)
			{
				BAF[IndxBAF1] <- rnorm(sum(IndxBAF1), mean=0.97, sd=0.01)
			}
		
			IndxBAF0 <- names(BAF) %in% names(Wave1PFB)[Wave1PFB < 0.05]
			if(sum(IndxBAF0) > 1)
			{
				BAF[IndxBAF0] <- rnorm(sum(IndxBAF0), mean=0.05, sd=0.01)
			}
	
			# Adding ramdom noise
			t  <- 1:length(X)
			ssp <- spectrum(X, plot=FALSE)  
			per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
			reslm <- lm(X ~ sin(2*pi/per*t)+cos(2*pi/per*t))		
			X <- X + (fitted(reslm)*20)
	
			# Adding bad SNPs (in general because of GC and LCR)
			TotalNumberofBadSNPs <- round(length(X)*BadSNPs[CHR])
			BadSNPsIndx <- sample(1:length(X), TotalNumberofBadSNPs)
			NoiseSNP <- sample(BadSNPIntensity, prob=BadSNPIntensityProb, 1)
			X[BadSNPsIndx] <- X[BadSNPsIndx] + rnorm(TotalNumberofBadSNPs, sd=(SD*2), mean=NoiseSNP)

			# BAF noise
			#BAF[BadSNPsIndx] <- BAF[BadSNPsIndx] + rnorm(TotalNumberofBadSNPs, sd=(SD), mean=0.1)
			
			# Add Telomere noise
			NTelomereSize <- sample(TelomereNoiseSize, 1)
			TeloEffect <- sample(TelomereNoiseEffect, 1) 
			X[1:NTelomereSize] <- X[1:NTelomereSize] + TeloEffect
			X[(length(X) - NTelomereSize):length(X)] <- X[(length(X) - NTelomereSize):length(X)] + TeloEffect
			
			# Adding CNVs		
			NumCNVs <- ((round(length(X)/1000))-2)
			DF <- sapply(1:NumCNVs, function(i) # Adding CNVs in the data.
			{
				CN <- sample(c(1,3), 1) # CNV Type
				PositionIndx <- as.numeric(i) * 1000
				#Size <- sample(CNVsSize, 1) # CNV Size
				Size <- CNVSizeFixed[i]
				
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
				
				# CNV position
				IndxV <- PositionIndx:(PositionIndx+Size)
				
				## Changing GLOBAL VARIABLES ##
				# LRR = X
				X[IndxV] <<- X[IndxV] + Impact
	
				# BAF, Change BAF but keep SNPs with low heterozygosity.
				NoChangeIndx <- c(which(BAF[IndxV] > 0.9), which(BAF[IndxV] < 0.1))
				NewIndx <- IndxV[(NoChangeIndx*-1)]
				BAF[NewIndx] <<- BAFCNV[NewIndx]
				## Changing GLOBAL VARIABLES ##
				BAF[BAF > 1] <<- 1
				BAF[BAF < 0] <<- 0
				
				df <- data.frame(Start=Position[PositionIndx], Stop=Position[(PositionIndx+Size)], StartIndx=PositionIndx, StopIndx=(PositionIndx+Size), NumSNPs=Size, Chr=CHR, CNVmean=Impact, CN=CN, sd=SD, ID=FileName, NoiseSNP=NoiseSNP, BadSNPs=TotalNumberofBadSNPs, NumCNVs=NumCNVs, stringsAsFactors=FALSE)
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
