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
##' @return Png files
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown. 
##' 

PlotHotspots <- function(CNVsDF, Hotspots, ListOfRawDataPath, Cores=1, Skip=10, Alpha=0.8, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18, OverlapMin=0.8, OverlapMax=1.2)
{
	suppressPackageStartupMessages(library(iPsychCNV))
	suppressPackageStartupMessages(library(mclust))
	suppressPackageStartupMessages(library(parallel))
	suppressPackageStartupMessages(library(biovizBase))
	suppressPackageStartupMessages(library(GenomeGraphs))
	suppressPackageStartupMessages(library(ggbio))
	suppressPackageStartupMessages(library(IRanges))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(GenomicRanges))
	suppressPackageStartupMessages(library(RColorBrewer))

	# loop over hotspots
	mclapply(1:nrow(Hotspots), mc.cores=Cores, mc.preschedule = FALSE, function(i) 
	{
		df <- SelectSamplesFromROI(DF=CNVsDF, roi=Hotspots[i,], OverlapMin=OverlapMin,  OverlapMax=OverlapMax)
		
		# Selecting only the best samples for hotspots plot
		if(nrow(df) > 20)
		{
			#df <- df[order(df$SDChr),]
			if(length(df$SDChr) > 0)
			{	
				del <- subset(df, CN == 1)
				dup <- subset(df, CN == 3)
				if(nrow(del) > 1){ del <- del[with(del, order(SDChr, SDCNV, CNVmean)), ] }
				if(nrow(dup) > 1){ dup <- dup[with(dup, order(SDChr, SDCNV, -CNVmean)), ] }
				
				if(nrow(del) > 10){ del <- del[1:10,] }
				if(nrow(dup) > 10){ dup <- dup[1:10,] }
				
				if(nrow(del) > 0 & nrow(dup) > 0)
				{
					df <- rbind(del,dup)
				}
				else if(nrow(del) == 0)
				{
					df <- dup
				}
				else
				{
					df <- del
				}
			}
			else
			{
				df <- df[1:20,]
			}
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
				CN <- X["CN"]
				if(Size < 100000){ Size <- Size * 2 }
	
				Start <- CNVstart - (Size)*(3/log10(NumSNP))^4
				Stop <-  CNVstop + (Size)*(3/log10(NumSNP))^4
		
				cat(ID, "\n")
				# Reading file
				File <- ListOfRawDataPath[grep(ID, ListOfRawDataPath)]
				Sample <- ReadSample(RawFile=File, skip=Skip, chr=Chr)
				Sample <- NormalizeData(Sample, ExpectedMean=0, penalty=penalty, Quantile=Quantile, QSpline=QSpline, sd=sd)
	
				# CNV region
				tmp <- subset(Sample, Position > Start & Position < Stop)
				tmp$CN <- CN
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
			
			#save(rect2, red2, file="Files.RData")
		
			p1 <- ggplot(red2, aes(Position, y = B.Allele.Freq)) + geom_point(aes(col=as.factor(CN)), size=0.5, alpha=Alpha) + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col="CNV region"), alpha=0.3, inherit.aes = FALSE) + theme(legend.title=element_blank()) + scale_color_manual(values = c("1" = Colors[1], "3"=Colors[3], "CNV region"=Colors[2]))

			# LogRRatio
			p2 <- ggplot(red2, aes(Position, y = Log.R.Ratio, col=as.factor(CN))) + geom_point(size=0.5, alpha=Alpha) + geom_line(aes(x=Position, y = Mean, col="Mean"), size = 0.5) + ylim(-1, 1) + theme(legend.title=element_blank()) + scale_color_manual(values = c("1" = Colors[1], "3"=Colors[3], "Mean"="black")) #  + scale_color_manual(values = c(Colors[1], "black"))
	
			Title <- paste("Hotspot: ", HotSpotID, sep="", collapse="")
			OutPlotfile <- paste("Hotspot", HotSpotID, "plot.png", sep=".", collapse="")
			Plot <- tracks(p3,p1, p2, main=Title, heights=c(3,5,5))
			ggsave(OutPlotfile, plot=Plot, height=5, width=10, dpi=200) 
		}
	})
}
