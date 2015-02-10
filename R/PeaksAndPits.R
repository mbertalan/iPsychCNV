##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PeaksAndPits
##' @return Data frame with QC variables.
##' @author Marcelo Bertalan
##' @export

PeaksAndPits <- function(DF, PlotName)
{
	library(pastecs)
	library(plyr)
	library(randtests)
	
	DF <- subset(DF, !Chr %in% c("MT", "X", "Y", "XY", "0"))
	tmp <- sapply(unique(DF$Chr), function(Z)
	{
		subCNV <- subset(DF, Chr %in% Z)
		
		subCNV$Log.R.Ratio[is.na(subCNV$Log.R.Ratio)] <- 0
		subCNV$B.Allele.Freq[is.na(subCNV$B.Allele.Freq)] <- 0
		BAF <- subCNV$B.Allele.Freq
		tmp <- smooth.spline(1:nrow(subCNV), subCNV$Log.R.Ratio)
		tp <- turnpoints(tmp$y)
		
		NewName <- paste(PlotName, "_chr", Z, ".png", sep="", collapse="")
			
		
		PointTest <- turning.point.test(subCNV$Log.R.Ratio)$p.value
		Peaks <- c(1,which(tp$peaks))
		Pits <- c(1,which(tp$pits))

		WaveLenPeaks <- sapply(2:length(Peaks), function(X){ Peaks[X] - Peaks[X-1] })
		WaveLenPits <- sapply(2:length(Pits), function(X){ Pits[X] - Pits[X-1] })

		
		WaveLength <- rep(0, length(tp$peaks))
		WaveLength[tp$peaks] <- WaveLenPeaks
		WaveLength[tp$pits] <- WaveLenPits
		Wave <- WaveLength[WaveLength > 0]

		Peaks <- c(1,which(tp$peaks))
		Pits <- c(1,which(tp$pits))
		Indx <- sort(c(Peaks,Pits))[-1]
		Points <- tp$points[Indx]
		Points2 <- Points + abs(min(Points))
		Amplitude <- sapply(1:(length(Points2)-1), function(X){ Points2[X] - Points2[X+1] })

		df <- data.frame(Position=tp$tppos, proba=tp$proba, WaveLength=Wave, Amplitude=Amplitude, stringsAsFactors=F)

		Log.R.Ratio.Median <- median(subCNV$Log.R.Ratio)
		Log.R.Ratio.Mean <- mean(subCNV$Log.R.Ratio)
		Log.R.Ratio.SD <- sd(subCNV$Log.R.Ratio)
		B.Allele.Freq.Median <- median(subCNV$B.Allele.Freq)
		B.Allele.Freq.Mean <- mean(subCNV$B.Allele.Freq)
		B.Allele.Freq.SD <- sd(subCNV$B.Allele.Freq)
		Wave.Length.Mean <- mean(df$WaveLength)
		Wave.Length.SD <- sd(df$WaveLength)
		Amplitude.Mean <- mean(df$Amplitude)
		Amplitude.SD <- sd(df$Amplitude)

		# BAF drift
		BAF_Res <- BAF_drift(BAF)
		BAF_Table <- table(BAF_Res$Class)
		BAF_Table <- (BAF_Table/sum(BAF_Table))*100
		BAF_Drift_Median <- median(BAF_Res$Drift)
		BAF_Drift_SD <- sd(BAF_Res$Drift)
		
		

		df2 <- data.frame(Log.R.Ratio.Median=Log.R.Ratio.Median, Log.R.Ratio.Mean=Log.R.Ratio.Mean, Log.R.Ratio.SD=Log.R.Ratio.SD, B.Allele.Freq.Median=B.Allele.Freq.Median,B.Allele.Freq.Mean=B.Allele.Freq.Mean, B.Allele.Freq.SD =B.Allele.Freq.SD , Wave.Length.Mean=Wave.Length.Mean, Wave.Length.SD=Wave.Length.SD, Amplitude.Mean =Amplitude.Mean , Amplitude.SD=Amplitude.SD, PointTest=PointTest, BAF_Drift_Median=BAF_Drift_Median, BAF_Drift_SD=BAF_Drift_SD, AAAA=BAF_Table["AAAA"], AAAB=BAF_Table["AAAB"], AAB=BAF_Table["AAB"],AB=BAF_Table["AB"], ABB=BAF_Table["ABB"], ABBB=BAF_Table["ABBB"], BBBB=BAF_Table["BBBB"], Chr=as.numeric(Z), ID=PlotName,stringsAsFactors=F)

		return(df2)
	})
	
	df <- apply(tmp, 2, function(X){ as.data.frame(X,  stringsAsFactors=F) })
	df <- ldply (df, data.frame)
	return(df)
}
