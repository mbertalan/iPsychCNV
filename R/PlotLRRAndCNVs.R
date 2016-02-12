##' PlotLRRAndCNVs: Plot  whole mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PlotLRRAndCNVs
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

PlotLRRAndCNVs <- function(CNV, Sample=MockData, CNVMean=0.3, Name="Test.png", Roi=RoiSingleMock)
{
	library(ggplot2)
	library(RColorBrewer)
	
	tmp <- Sample
	
	if(length(CNV$Source) == 0)
	{
		CNV$Source <- "iPsychCNV" 
	}
	
	if(length(CNV$CNVmean) == 0)
	{
		if(length(CNV$CNVMean) == 0)
		{
			CNV$CNVmean <- rep(CNVMean, nrow(CNV))
			CNV$CNVmean[CNV$CN == 1] <- CNV$CNVmean[CNV$CN == 1] * -1
			CNV$CNVmean[CNV$CN == 0] <- Roi$CNVmean[CNV$CN == 0] * -1
		}
		else
		{
			CNV$CNVmean <- CNV$CNVMean
		}
	}
	
	if(length(Roi$CNVmean) == 0)
	{
		Roi$CNVmean <- rep(CNVMean, nrow(Roi))
		Roi$CNVmean[Roi$CN == 1] <- Roi$CNVmean[Roi$CN == 1] * -1
		Roi$CNVmean[Roi$CN == 0] <- Roi$CNVmean[Roi$CN == 0] * -1
	}
	Roi$YMin <- rep(0, length(Roi$CNVmean))
	Roi$YMax <- rep(1, length(Roi$CNVmean))
	
	Mean <- SlideWindowMean(tmp$Log.R.Ratio, 35)
	tmp$Mean <- Mean


	Colors <- brewer.pal("Set1", n=9)
	Colors2 = brewer.pal("Dark2", n=7)

	# LRR
	p1 <- ggplot(tmp, aes(x=Position, y=Log.R.Ratio)) + geom_point(aes(x=Position, y=Log.R.Ratio), alpha=0.05, size=1) 
	p1 <- p1 + geom_rect(data=Roi, aes(xmin =Start, xmax=Stop, ymin=(CNVmean-0.1), ymax=(CNVmean+0.1)), colour=Colors2[3], fill=NA, alpha=0.3, inherit.aes = FALSE, size=0.5) + theme(legend.title=element_blank())
	p1 <- p1 +  geom_segment(data=CNV, aes(x = Start, y = CNVmean, xend = Stop, yend = CNVmean, colour=as.factor(CN)), size=4) + scale_colour_manual(values = c("1" = Colors[1], "2"=Colors[9], "3" = Colors[2], "4" = Colors[3], "0"=Colors[4], "5"=Colors[5], "B.Allele.Freq"=Colors[2], "CNV region"=Colors[3], "CNV predicted"=Colors[4], "Mean"="black"))	
	p1 <- p1 + geom_line(data=tmp, aes(x=Position, y = Mean), size = 0.1, alpha=0.2) 
	
	# BAF
	p2 <- ggplot(tmp, aes(x=Position, y=B.Allele.Freq)) + geom_point(aes(x=Position, y = B.Allele.Freq, col="B.Allele.Freq"), alpha=0.3, size=1) 
	p2 <- p2 + geom_rect(data=Roi, aes(xmin=Start, xmax=Stop, ymin=YMin, ymax=YMax, col="CNV region"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) 
	
	CNV <- subset(CNV, CN != 2)
	if(nrow(CNV) > 0)
	{
		retCNV <- data.frame(Source=CNV$Source, Start=CNV$Start, Stop=CNV$Stop, ymin=rep(0.4, length(CNV$Stop)), ymax=rep(0.6, length(CNV$Stop)))
		p2 <- p2 + geom_rect(data=retCNV, aes(xmin=Start, xmax=Stop, ymin=ymin, ymax=ymax, col="CNV predicted"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) 
	}

	p2 <- p2 + scale_colour_manual(values = c("0"=Colors[4], "1" = Colors[1], "2"=Colors[9], "3" = Colors[2], "4" = Colors[3], "5" = Colors[5] , "B.Allele.Freq"=Colors[2], "CNV region"=Colors[3], "CNV predicted"=Colors[4], "Mean"="black"))		
	library(ggbio)
	Plot <- tracks(p1, p2, heights=c(4,4))
	ggsave(Plot, file=Name, width=10, height=6, dpi=300)
	#
}
