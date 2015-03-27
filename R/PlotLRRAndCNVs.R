##' PlotLRRAndCNVs: Plot  whole mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PlotLRRAndCNVs
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

PlotLRRAndCNVs <- function(PennCNV, tmp=MockData, CNVMean, Name="Test.png", Roi=RoiSingleMock)
{
	library(ggplot2)
	library(ggbio)

	Roi$CNVmean <- rep(CNVMean, 24)
	Roi$CNVmean[Roi$CN == 1] <- Roi$CNVmean[Roi$CN == 1] * -1
	Roi$CNVmean[Roi$CN == 0] <- Roi$CNVmean[Roi$CN == 0] * -1
	Roi$YMin <- rep(0, length(Roi$CNVmean))
	Roi$YMax <- rep(1, length(Roi$CNVmean))
	
	if(!length(PennCNV$CNVmean) > 0)
	{
		PennCNV$CNVmean <- rep(CNVMean, nrow(PennCNV))
		PennCNV$CNVmean[PennCNV$CN == 1] <- PennCNV$CNVmean[PennCNV$CN == 1] * -1
		PennCNV$CNVmean[PennCNV$CN == 0] <- PennCNV$CNVmean[PennCNV$CN == 0] * -1
	}
	
	Mean <- SlideWindowMean(tmp$Log.R.Ratio, 35)
	tmp$Mean <- Mean


	Colors <- brewer.pal("Set1", n=9)
	Colors2 = brewer.pal("Dark2", n=7)

	# LRR
	p1 <- ggplot(tmp, aes(Position, Log.R.Ratio))
	p1 <- p1 + geom_point(alpha=0.05, size=1) 
	p1 <- p1 + geom_rect(data=Roi, aes(xmin =StartPos, xmax=StopPos, ymin=(CNVmean-0.1), ymax=(CNVmean+0.1)), colour=Colors2[3], fill=NA, alpha=0.3, inherit.aes = FALSE, size=0.5) + theme(legend.title=element_blank())
	p1 <- p1 +  geom_segment(data=PennCNV, aes(x = Start, y = CNVmean, xend = Stop, yend = CNVmean, colour=as.factor(CN)), size=4) + scale_colour_manual(values = c("1" = Colors[1], "2"=Colors[9], "3" = Colors[2], "4" = Colors[3], "0"=Colors[4], "B.Allele.Freq"=Colors[2], "CNV region"=Colors[3], "CNV predicted"=Colors[4], "Mean"="black"))	
	p1 <- p1 + geom_line(data=tmp, aes(x=Position, y = Mean), size = 0.1, alpha=0.4) 
	
	# BAF
	p2 <- ggplot(tmp, aes(Position, y = B.Allele.Freq)) 
	p2 <- p2 + geom_point(aes(col="B.Allele.Freq"), alpha=0.6, size=1)  
	p2 <- p2 + geom_rect(data=Roi, aes(xmin=StartPos, xmax=StopPos, ymin=YMin, ymax=YMax, col="CNV region"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) 
	
	PennCNV <- subset(PennCNV, CN != 2)
	if(nrow(PennCNV) > 0)
	{
		retPennCNV <- data.frame(Start=PennCNV$Start, Stop=PennCNV$Stop, ymin=rep(0.4, length(PennCNV$Stop)), ymax=rep(0.6, length(PennCNV$Stop)))
		p2 <- p2 + geom_rect(data=retPennCNV, aes(xmin=Start, xmax=Stop, ymin=ymin, ymax=ymax, col="CNV predicted"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) 
	}

	p2 <- p2 + scale_colour_manual(values = c("1" = Colors[1], "2"=Colors[9], "3" = Colors[2], "4" = Colors[3], "0"=Colors[4], "B.Allele.Freq"=Colors[2], "CNV region"=Colors[3], "CNV predicted"=Colors[4], "Mean"="black"))		

	Plot <- tracks(p1, p2, heights=c(4,4))
	ggsave(Plot, file=Name, width=16, height=6, dpi=300)
	#
}
