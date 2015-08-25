##' PlotHotspots
##'
##' @title PlotHotspots
##' @return png files
##' @author Marcelo Bertalan
##' @export

PlotHotspots <- function(CNVsDF, Hotspots, AllHeadFiles, Cores=1, Skip=10)
{
	library(iPsychCNV)
	library(mclust)
	library(parallel)
	library(biovizBase)
	library(GenomeGraphs)
	library(ggbio)
	library(IRanges)
	library(ggplot2)
	library(GenomicRanges)
	library(RColorBrewer)

	# loop over hotspots
	mclapply(1:nrow(Hotspots), mc.cores=Cores, mc.preschedule = FALSE, function(i) 
	{
		df <- SelectSamplesFromROI(DF=CNVsDF, roi=Hotspots[i,], OverlapMin=0.9,  OverlapMax=1.1)
		
		# Selecting only the best samples for hotspots plot
		if(nrow(df) > 10)
		{
			df <- df[order(df$SDChr),]
			df <- df[1:10,]
		}

		if(nrow(df) > 0)
		{
			HotSpotID <- paste(Hotspots[i,"Chr"],Hotspots[i,"Start"], Hotspots[i,"Stop"], sep=".", collapse="")
			Start <- Hotspots[i,"Start"]
			Stop <- Hotspots[i,"Stop"]	

			# collecting data
			Data <- apply(df, 1, function(X)
			{
				Chr <- X["Chr"]
				CNVstart <- as.numeric(X["Start"])
				CNVstop <- as.numeric(X["Stop"])
				ID <- X["ID"]
				NumSNP <- as.numeric(X["NumSNPs"])
				Size <- as.numeric(X["Length"])
	
				Start <- CNVstart - (Size)*(3/log10(NumSNP))^3
				Stop <-  CNVstop + (Size)*(3/log10(NumSNP))^3
		
				cat(ID, "\n")
				# Reading file
				File <- AllHeadFiles[grep(ID, AllHeadFiles)]
				Sample <- ReadSample(RawFile=File, skip=Skip, chr=Chr)
	
				# CNV region
				tmp <- subset(Sample, Position > Start & Position < Stop)
				return(tmp)
			})
			red <- MatrixOrList2df(Data)
			red2 <- red[with(red, order(Position)),]	
			Mean <- SlideWindowMean(red2$Log.R.Ratio, 35)
			red2$Mean <- Mean	
			chr <- unique(red2$Chr)[1]
				
			# Ideogram
			data(hg19IdeogramCyto,package="biovizBase")
			CHR <- paste("chr", chr, collapse="", sep="")
			p3 <- plotIdeogram(hg19IdeogramCyto,CHR,cytoband=TRUE,xlabel=TRUE, aspect.ratio = 1/85, alpha = 0.3) + xlim(GRanges(CHR,IRanges(Start,Stop)))
	
			# Colors
			Colors = brewer.pal(7,"Set1")
	
			# B.Allele
			rect2 <- data.frame (xmin=Start, xmax=Stop, ymin=0, ymax=1) # CNV position
		
			p1 <- ggplot(red2, aes(Position, y = B.Allele.Freq)) + geom_point(aes(col="B.Allele.Freq"), size=1, alpha = 0.3) + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col="CNV region"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[2:4]))

			# LogRRatio
			p2 <- ggplot(red2, aes(Position, y = Log.R.Ratio, col="Log.R.Ratio")) + geom_point(alpha = 0.3, size=1) + geom_line(aes(x=Position, y = Mean, col="Mean"), size = 0.5) + ylim(-1, 1) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[1], "black"))
	
			Title <- paste("Hotspot: ", HotSpotID, sep="", collapse="")
			OutPlotfile <- paste("Hotspot", HotSpotID, "plot.png", sep=".", collapse="")
			Plot <- tracks(p3,p1, p2, main=Title, heights=c(3,5,5))
			ggsave(OutPlotfile, plot=Plot, height=5, width=10, dpi=100) 
		}
	})
}
