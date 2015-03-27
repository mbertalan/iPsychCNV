PlotCNVsFromDataFrame <- function(DF, Cores) # Path from cluster  
{
	library(mclust)
	library(parallel)

	DF$UniqueID <- 1:nrow(DF)
	matrixList <- lapply(1:NROW(as.matrix(DF)), function(i) DF[i,,drop=FALSE]) 
	mclapply(matrixList, mc.cores=Cores, function(X) 
	{
		
		NAMES <- names(X)
		X2 <- unlist(X)
		names(X2) <- NAMES
		X <- X2
		tmpSize <- nrow(X)
		tmpLength <- length(X)
		
		chr <- as.numeric(X["Chr"])
		ID <- X["ID"]
		UniqueID <- X["UniqueID"]

		cat(ID, "\n")
	
		CNVstart <- as.numeric(X["Start"])
		CNVstop <- as.numeric(X["Stop"])
		Size <- as.numeric(X["Length"])
		#Source <- X["Source"]
		CN <- X["Class"]
		SDCNV <- round(as.numeric(X["SDCNV"]), digits=2)
		NumSNP <- as.numeric(X["NumSNPs"])

		#if("PlotStart" %in% names(X))
		#{
		#	Start <- as.numeric(X["PlotStart"])
		#	Stop <- as.numeric(X["PlotStop"])
		#}
		#else
		#{
			Start <- CNVstart - Size*(3/log10(NumSNP))^3
			Stop <-  CNVstop + Size*(3/log10(NumSNP))^3
		#}
		Conf <- 100
		NewName <- paste(UniqueID, CNVstart, CNVstop, chr, ID, sep="_", collapse="")
		OutPlotfile <- paste(NewName, "plot.png", sep=".", collapse="")
 		#cat(ID, CN, "\n")
		# Version 3
		# Include Overlap cutoff
		# Relative size to CNV
		# Start and stop = Region to plot
		# CNVstart, CNVstop, = CNV region.
	
		library(biovizBase)
		library(GenomeGraphs)
		library(ggbio)
		library(IRanges)
		library(ggplot2)
		library(GenomicRanges)
		library(RColorBrewer)

		# Reading sample file
		source("~/CNV/Scripts/iPsychCNV/CNV_Functions.R")
		PathRawData <- GetSamplePath(ID)
		RawFile <- paste(PathRawData, ID, sep="", collapse="")
		cat(RawFile,"\n")
		sample<-read.table(RawFile,as.is=T,head=T,sep="\t", skip=10)
		red<-subset(sample,Chr==chr)
		save(red, file="red.RData")	
		red <-subset(red, Position > Start & Position < Stop)	
		red$Log.R.Ratio[is.na(red$Log.R.Ratio)] <- 0

		red2 <- red[with(red, order(Position)),]
		save(red2, file="red2.RData")	
		red2$Log.R.Ratio[is.na(red2$Log.R.Ratio)] <- 0
		source("~/CNV/Scripts/iPsychCNV/SlidingWindowMean.R")
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
