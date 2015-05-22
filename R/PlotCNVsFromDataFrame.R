##' Plot CNVs from data frame 
##'
##' @title PlotCNVsFromDataFrame
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

PlotCNVsFromDataFrame <- function(DF, PathRawData=".", Cores=1, Skip=0, PlotPosition=1, Pattern="*",recursive=TRUE) # Path from cluster  
{
	library(mclust)
	library(parallel)
	library(biovizBase)
	library(GenomeGraphs)
	library(ggbio)
	library(IRanges)
	library(ggplot2)
	library(GenomicRanges)
	library(RColorBrewer)

	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
	DF$UniqueID <- 1:nrow(DF)
	#matrixList <- lapply(1:NROW(as.matrix(DF)), function(i) DF[i,,drop=FALSE]) 
	#mclapply(matrixList, mc.cores=Cores, function(X) 
	mclapply(DF$UniqueID, mc.cores=Cores, function(UID) 
	{
		X <- subset(DF, UniqueID %in% UID)
	
		chr <- X$Chr
		ID <- X$ID
		UniqueID <- X$UniqueID

		cat(ID, "\n")
	
		CNVstart <- as.numeric(X$Start)
		CNVstop <- as.numeric(X$Stop)
		Size <- as.numeric(X$Length)
		CN <- X$CN
		SDCNV <- round(as.numeric(X$SDCNV), digits=2)
		NumSNP <- as.numeric(X$NumSNPs)

		Start <- CNVstart - (Size*PlotPosition)*(3/log10(NumSNP))^3
		Stop <-  CNVstop + (Size*PlotPosition)*(3/log10(NumSNP))^3
	
		NewName <- paste(UniqueID, CNVstart, CNVstop, chr, ID, sep="_", collapse="")
		OutPlotfile <- paste(NewName, "plot.png", sep=".", collapse="")

		# Reading sample file
		#RawFile <- paste(PathRawData, ID, sep="", collapse="")
		RawFile <- Files[grep(ID, Files)]
		cat(RawFile,"\n")
		sample <- ReadSample(RawFile, skip=Skip)
		red<-subset(sample,Chr==chr)	
		red <-subset(red, Position > Start & Position < Stop)	
		red2 <- red[with(red, order(Position)),]
		
	
		Mean <- SlideWindowMean(red2$Log.R.Ratio, 35)
		red2$Mean <- Mean
		

		# Ideogram
		data(hg19IdeogramCyto,package="biovizBase")
		CHR <- paste("chr", chr, collapse="", sep="")
		p3 <- plotIdeogram(hg19IdeogramCyto,CHR,cytoband=TRUE,xlabel=TRUE, aspect.ratio = 1/85, alpha = 0.3) + xlim(GRanges(CHR,IRanges(Start,Stop)))
	
		# Colors
		Colors = brewer.pal(7,"Set1")

		# B.Allele
		rect2 <- data.frame (xmin=CNVstart, xmax=CNVstop, ymin=0, ymax=1) # CNV position
	
		p1 <- ggplot(red2, aes(Position, y = B.Allele.Freq)) + geom_point(aes(col="B.Allele.Freq"), size=1) + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col="CNV region"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[2:4]))

		# LogRRatio
		p2 <- ggplot(red2, aes(Position, y = Log.R.Ratio, col="Log.R.Ratio")) + geom_point(alpha = 0.6, size=1) + geom_line(aes(x=Position, y = Mean, col="Mean"), size = 0.5) + ylim(-1, 1) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[1], "black"))

		Title <- paste("CN: ", CN, ", SDCNV: ", SDCNV, ", NumSNPs: ", NumSNP, ", Sample: ", ID, sep="", collapse="")
		Plot <- tracks(p3,p1, p2, main=Title, heights=c(3,5,5))
		ggsave(OutPlotfile, plot=Plot, height=5, width=10, dpi=100) 
	})
}
