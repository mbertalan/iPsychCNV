##' PlotCNVsFromDataFrame: Function to plot Log R Ratio (LRR) and B Allele Frequency (BAF) of CNVs from a data frame (DF). 
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##'
##' @title PlotCNVsFromDataFrame
##' @param DF: Unknown, default = Unknown. 
##' @param PathRawData: Path for the raw data. 
##' @param Cores: Unknown, default = 1.
##' @param Skip: Unknown, default = 0.
##' @param PlotPosition: Unknown, default = 1.
##' @param Pattern: Unknown, default = "*".
##' @param Recursive: Unknown, default = TRUE.
##' @param Dpi: Dots per inch, default = 300. 
##' @param Files: Dots per inch, default = NA.
##' @param SNPList: Getting chromosome (Chr) and position from another source than the RawFile - input should be the full path of the SNPList with columns: Name, Chr, and Position. Any positions from the RawFile will be erased. A PFB-column is also allowed but will be overwritten by the PFB-parameter or exchanged with 0.5, default = NULL.   
##' @param Key: Exchange the ID printed on the plot and in the name of file with a deidentified ID - requires that the DF contains a column called ID_deidentified, default = NA. 
##' @param OutFolder: Path for saving outputfiles, default is the current folder.
##' @return One BAF- and LRR-plot for each CNV.
##' @author Marcelo Bertalan, Ida Elken Sønderby, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples 
##'PlotCNVsFromDataFrame(DF=CNVs.Good[1:4,], PathRawData=".", Cores=1, Skip=0, Pattern="^MockSamples*", key=NA, OutFolder=".")


PlotCNVsFromDataFrame <- function(DF, PathRawData=".", Cores=1, Skip=0, PlotPosition=1, Pattern="*",recursive=TRUE, dpi=300, Files=NA, SNPList=NULL, key=NA, OutFolder=".") # Path from cluster  ByFolder=FALSE, 
{
  library(ggplot2)
  library(ggbio) # For some reason ggplot2 2.0.2 is not working, probably conflict with other packages. Version 1.0.1 works.
  library(parallel)
  library(biovizBase)
  library(RColorBrewer)
  library(GenomicRanges)
  
  LocalFolder <- PathRawData
  if(is.na(Files))
  {
    Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
  }
  
  DF$UniqueID <- 1:nrow(DF)
  
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

    ## Naming output-file    
 
    # based on key or not
    if (!is.na(key)) { # if want a different ID from the genetic ID in the plot
      ID_deidentified <- X$ID_deidentified  # Added this to get ID_deidentified
      NewName <- paste(ID_deidentified, "_", UniqueID, "_chr", chr, ":", CNVstart, "-", CNVstop, sep="", collapse="")
    } else {
      NewName <- paste(ID, "_", UniqueID, "_chr", chr, ":", CNVstart, "-", CNVstop, sep="", collapse="")
    }

    # based on OutFolder or not  
      if(OutFolder!=".") {
      OutPlotfile <- paste(OutFolder, NewName, "_plot.png", sep="", collapse="")
      print(OutPlotfile)
      }
    else {
      OutPlotfile <- paste(NewName, "plot.png", sep="_", collapse="")
      print(OutPlotfile)
      }

    # Reading sample file
    #RawFile <- paste(PathRawData, ID, sep="", collapse="")
    RawFile <- Files[grep(ID, Files)]
    cat("File: ", RawFile,"\n")
    
    sample <- ReadSample(RawFile, skip=Skip, SNPList=SNPList)
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
    
    p1 <- ggplot(red2, aes(Position, y = B.Allele.Freq)) + geom_point(aes(col="B.Allele.Freq"), size=0.5) + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col="CNV region"), alpha=0.2, inherit.aes = FALSE) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[2:4]))
    
    # LogRRatio
    p2 <- ggplot(red2, aes(Position, y = Log.R.Ratio, col="Log.R.Ratio")) + geom_point(alpha = 0.6, size=0.5) + geom_line(aes(x=Position, y = Mean, col="Mean"), size = 0.5) + ylim(-1, 1) + theme(legend.title=element_blank()) + scale_color_manual(values = c(Colors[1], "black"))
    
    # Title printed for plot
    if (!is.na(key)) { # if want a different ID from the genetic ID in the plot
      Title <- paste("CN: ", CN, ", SDCNV: ", SDCNV, ", NumSNPs: ", NumSNP, ", Sample: ", ID_deidentified, sep="", collapse="")
    } else {
      Title <- paste("CN: ", CN, ", SDCNV: ", SDCNV, ", NumSNPs: ", NumSNP, ", Sample: ", ID, sep="", collapse="")
    }
        Plot <- tracks(p3,p1, p2, main=Title, heights=c(3,5,5))
    ggsave(OutPlotfile, plot=Plot, height=5, width=10, dpi=dpi) 
  })
}
