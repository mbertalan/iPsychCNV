##' PlotHotspots: Plot hot spot for CNVs.
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PlotHotspots
##' @param CNVsDF: Data frame with CNVs, default = Unknown.
##' @param Hotspots: Unknown, default = Unknown.
##' @param ListOfRawDataPath: Full path for each sample. You can get it using the command: Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive). Than use ListOfRawDataPath=Files.
##' @param Cores: The number of cores used, default = 1.
##' @param Skip: Integer, the number of lines of the data file to be skipped before beginning to read the data, default = 10.
##' @param Penalty: The coefficient of the penalty for degrees of freedom, default = 60.
##' @param Quantile: Logical, if quantile normalization should be applied or not, default = FALSE.
##' @param QSpline: Logical, if a cubic smoothing spline should be used to normalize the data, default = FALSE.
##' @param Sd: Numeric, Log R Ratio standard deviation for the quantile normalization, default = 0.18.
##' @param SNPList: Getting chromosome (Chr) and position from another source than the RawFile - input should be the full path of the SNPList with columns: Name, Chr, and Position. Any positions from the RawFile will be erased. A PFB-column is also allowed but will be overwritten by the PFB-parameter or exchanged with 0.5, default = NULL.
##' @return Png files
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

PlotHotspots <- function(CNVsDF, Hotspots, ListOfRawDataPath, Cores=1, Skip=10, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18, OverlapMin=0.8, OverlapMax=1.2)
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
		df <- SelectSamplesFromROI(DF=CNVsDF, roi=Hotspots[i,], OverlapMin=OverlapMin,  OverlapMax=OverlapMax)

		# Selecting only the best samples for hotspots plot
		if(nrow(df) > 10)
		{
			#df <- df[order(df$SDChr),]
			if(length(df$SDChr) > 0)
			{
				df <- df[with(df, order(SDChr, !CN)), ]
			}
			df <- df[1:10,]
		}

		if(nrow(df) > 0)
		{
			HotSpotID <- paste(Hotspots[i,"Chr"],Hotspots[i,"Start"], Hotspots[i,"Stop"], sep=".", collapse="", SNPList = NULL)
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
				if(Size < 100000){ Size <- Size * 2 }

				Start <- CNVstart - (Size)*(3/log10(NumSNP))^4
				Stop <-  CNVstop + (Size)*(3/log10(NumSNP))^4

				cat(ID, "\n")
				# Reading file
				File <- ListOfRawDataPath[grep(ID, ListOfRawDataPath)]
				Sample <- ReadSample(RawFile=File, skip=Skip, chr=Chr, SNPList = SNPList)
				Sample <- NormalizeData(Sample, ExpectedMean=0, penalty=penalty, Quantile=Quantile, QSpline=QSpline, sd=sd)

				# CNV region
				tmp <- subset(Sample, Position > Start & Position < Stop)
				return(tmp)
			})
			red <- MatrixOrList2df(Data)
			red2 <- red[with(red, order(Position)),]
			WindowSize <- round(nrow(red2)/100)
			if(WindowSize < 5){ WindowSize <- 5 }
			if(WindowSize > 35){ WindowSize <- 35 }
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

			p1 <- ggplot(red2, aes(Position, y = B.Allele.Freq)) + geom_point(aes(col="B.Allele.Freq"), size=1) + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col="CNV region"), alpha=0.3, inherit.aes = FALSE) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[2:4]))

			# LogRRatio
			p2 <- ggplot(red2, aes(Position, y = Log.R.Ratio, col="Log.R.Ratio")) + geom_point(size=1) + geom_line(aes(x=Position, y = Mean, col="Mean"), size = 0.5) + ylim(-1, 1) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[1], "black"))

			Title <- paste("Hotspot: ", HotSpotID, sep="", collapse="")
			OutPlotfile <- paste("Hotspot", HotSpotID, "plot.png", sep=".", collapse="")
			Plot <- tracks(p3,p1, p2, main=Title, heights=c(3,5,5))
			ggsave(OutPlotfile, plot=Plot, height=5, width=10, dpi=200)
		}
	})
}
