##' PlotAllCNVs: Plot all Copy Number Variations (CNVs) in one plot. 
##'
##' Specifically designed to handle noisy data from amplified DNA on Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PlotAllCNVs
##' @param Df: Data frame with predicted CNVs for each sample, default = CNV.Res.
##' @param Name:Unknown, default = Unknown.
##' @param NCOL: Number of columns, default = 2.
##' @param Roi: Regions of interest or hotspots, default = roi.
##' @param Width: Unknown, defualt = 16.
##' @param Height: Unknown, default = 30.
##' @param Hg: Human genome version, default = hg19.
##' @param Chr: Chromosome, default = NA.
##' @param Start: The start position of each CNV, default = NA.
##' @param Stop: The stop position of each CNV, default = NA.
##' @param Order: Should the function order the CNVs ?, default = TRUE.
##' @return Return data in a data frame.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export 
##' @examples Unknown.
##' 

PlotAllCNVs <- function(df=CNV.Res, Name="CNV.Res.Test.png", NCOL=2, Roi=roi, width=16, height=30, hg="hg19", chr=NA, start=NA, stop=NA, Order=TRUE) # TimesLength
{	
	library(scales)
	library(ggplot2)
	library(RColorBrewer)
	library(grid)
	
	# Select which human build it will plot
	CytoBands <- subset(CytoBands2, Human %in% hg & Chr %in% c(1:22))	

	if(length(df$CN) == 0){ stop("Missing CN\n") }
	
	if(length(df$CN) > 0){ df$Class <- df$CN }
	if(length(df$File) > 0){ df$ID <- df$File }

	tmp <- df[, c("Start","Stop","Chr","Length","ID", "Class", "CN")]
	tmp <- subset(tmp, Chr %in% c(1:22))
	tmp3 <- sapply(unique(tmp$Chr), function(X)
	{
		tmp2 <- subset(tmp, Chr %in% X)
		# Should I order your CNVs ?
		if(Order){ tmp2 <- tmp2[with(tmp2, order(Start, !Length, CN)), ]  })
		Samples <- 1:length(unique(tmp2$ID))
		names(Samples) <- unique(tmp2$ID)
		Indx <- sapply(tmp2$ID, function(Y){ Samples[as.character(Y)] })
		tmp2$Indx <- Indx
		return(tmp2)
	})	
	tmp2 <- MatrixOrList2df(tmp3)
	
	# Join df + roi + cytobands 
	tmp2 <- tmp2[, c("Start","Stop","Chr","Length","ID", "Class", "Indx")]
	roi <- Roi[,c("Start","Stop","Chr","Length","ID", "Class")]
	CytoBands <- CytoBands[,c("Start","Stop","Chr","Length","ID", "Class")]

	CytoBands$Indx <- rep(-1, nrow(CytoBands))
	roi$Indx <- rep(-2, nrow(roi))
	
	tmp2 <- rbind(tmp2, CytoBands, roi)

	# Subset if specific region is selected.
	if(!is.na(chr))
	{ 
		tmp2 <- subset(tmp2, Chr %in% chr & Start < stop & Stop > start)
	}	

	# Creating titles and check CNV count.
	tmp3 <- subset(tmp2, Indx > 0)
	Info <- tapply(tmp3$ID, as.factor(tmp3$Chr), function(X){ c(length(unique(X)), length(X)) }) 
	df <- do.call(rbind.data.frame, Info)
	names(df) <- c("Samples", "CNVs")
	df$Samples <- as.numeric(df$Samples) 
	df$CNVs <- as.numeric(df$CNVs)
	df$Chr <- rownames(df)
	
	Title <- apply(df, 1, function(X){ paste("chr:", X["Chr"], ", Samples:", X["Samples"], ", CNVs:", X["CNVs"], sep="", collapse="") })
	Titles <- sapply(tmp2$Chr, function(X){ if(is.na(Title[X])){ paste("chr:", X, ", Samples: 0, CNVs: 0", sep="", collapse="") }else{ Title[X]  } })
	tmp2$Title <- Titles
	
	TitleIndx <- sapply(unique(tmp2$Title), function(X){ as.numeric(unlist(strsplit(unlist(strsplit(unique(X)[1], ","))[1], "chr:"))[2]) })
	tmp2$Titles <-  factor(tmp2$Title, levels= names(sort(TitleIndx)))

	tmp2$Start <- as.numeric(tmp2$Start)
	tmp2$Stop <- as.numeric(tmp2$Stop)

	tmp2ROI <- subset(tmp2, Class %in% "ROI")
	# Change name from start and stop so it is different from the two data.frames. This avoid problems with ggplot2 aes.
	
	Colors <- brewer.pal("Set1", n=9) # + scale_color_brewer(palette="Set1")
	b <- ggplot() + geom_segment(data=tmp2, aes(x = Start, y = Indx, xend = Stop, yend = Indx, colour=as.factor(Class)))# aes(x=Start, y=Indx)
	b <- b + scale_colour_manual(values = c("ROI" = Colors[6],"q" = Colors[5],"p" = Colors[4], "1" = Colors[1], "3" = Colors[2], "4" = Colors[3], "0"=Colors[7], "2"=Colors[9], "5" = Colors[5]))
	b <- b + geom_vline(aes(xintercept = Start), data=tmp2ROI, alpha=0.2) 
	b <- b + geom_vline(aes(xintercept = Stop), data=tmp2ROI, alpha=0.2)
	
	if(is.na(start))
	{
		b <- b + facet_wrap(~ Titles, scales = "free", ncol = NCOL) 
	}else{
		b <- b + xlim(start, stop)
	}
	
	ggsave(b, file=Name, width=width, height=height, dpi=300)
}
