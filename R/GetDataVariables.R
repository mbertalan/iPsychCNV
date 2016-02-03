##' GetDataVariables
##'
##' @title CompressCNVsHighFreq
##' @return Hotspots
##' @author Marcelo Bertalan
##' @export

GetDataVariables <- function(Data)
{
	library(randtests)
	library(moments)
	library(stats)
	suppressPackageStartupMessages(library(pastecs))
	
	# Score <- ((abs(AMP.df40.good$CNVmeanNorm)*10)*(BothDiff*10))/(AMP.df40.good$AB+1)
	ANS <- 1; AGO <- 1; BOX <- 1; COX <- 1; 
	
	Window <- round(length(Data$CNV)/50)
	if(Window < 3){ Window <- 3 }
	MeansCNV <- SlideWindowMean(Data$LRR, Window)
	SDMeansCNV <- sd(MeansCNV)
	CNVmeanOrig <- mean(Data$LRR) # Original, not norm.
	CNVmean <- mean(Data$CNV) 
	HighMean <- mean(Data$High)
	LowMean <- mean(Data$Low)
	VarianceCNV <- var(Data$LRR)
	SDCNV <- sd(Data$LRR)
	SDHigh <- sd(Data$High)
	SDLow <- sd(Data$Low)
	DiffHigh <- abs(CNVmean - HighMean)
	DiffLow <- abs(CNVmean - LowMean)
	CNVmeanByRef <- CNVmean - ((HighMean+LowMean)/2)
	
	# Wave
	#tp <- turnpoints(Data$LRR)
	#Peaks <- c(1,which(tp$peaks))
	#Pits <- c(1,which(tp$pits))

	#WaveLenPeaks <- sapply(2:length(Peaks), function(X){ Peaks[X] - Peaks[X-1] })
	#WaveLenPits <- sapply(2:length(Pits), function(X){ Pits[X] - Pits[X-1] })
		
	#WaveLength <- rep(0, length(tp$peaks))
	#WaveLength[tp$peaks] <- WaveLenPeaks
	#WaveLength[tp$pits] <- WaveLenPits
	#Wave <- WaveLength[WaveLength > 0]

	#Indx <- sort(c(Peaks,Pits))[-1]
	#Points <- tp$points[Indx]
	#Points2 <- Points + abs(min(Points))
	#Amplitude <- sapply(1:(length(Points2)-1), function(X){ Points2[X] - Points2[X+1] })
	#Wave.Length.Mean <- mean(Wave)
	#Wave.Length.SD <- sd(Wave)
	#Amplitude.Mean <- mean(Amplitude)
	#Amplitude.SD <- sd(Amplitude)

	#res2 <- data.frame(CNVmean=CNVmean, HighMean=HighMean, LowMean=LowMean, CNVmeanOrig=CNVmeanOrig, SDCNV=SDCNV, SDMeansCNV=SDMeansCNV, SDHigh=SDHigh, SDLow=SDLow,VarianceCNV=VarianceCNV, DiffHigh=DiffHigh, DiffLow=DiffLow,  Wave.Length.Mean=Wave.Length.Mean, Wave.Length.SD=Wave.Length.SD, Amplitude.Mean=Amplitude.Mean, Amplitude.SD=Amplitude.SD)
	res2 <- data.frame(CNVmean=CNVmean, CNVmeanOrig=CNVmeanOrig, CNVmeanByRef=CNVmeanByRef, UpMean=HighMean, DownMean=LowMean, SDChr=SDCNV) # SDMeansCNV=SDMeansCNV, SDHigh=SDHigh, SDLow=SDLow,VarianceCNV=VarianceCNV, DiffHigh=DiffHigh, DiffLow=DiffLow, CNVmeanOrig=CNVmeanOrig
	return(res2)	
	
}
