##' PlotAllInOnePlotggplot: Plot all CNVs in one plot. 
##'
##' @title PlotAllinOnePlotggplot
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

PlotAllInOnePlotggplot <- function(tmp, Name="Test.png", NCOL=2, roi, width=16, height=30) # TimesLength
{	
	library(scales)
	library(ggplot2)
	library(RColorBrewer)

	if(length(tmp$Class) > 0){ tmp$CN <- tmp$Class }
	if(length(tmp$File) > 0){ tmp$ID <- tmp$File }

	tmp <- tmp[, c("Start","Stop","Chr","Length","ID", "Class")]
	tmp <- subset(tmp, Chr %in% c(1:22))
	roi <- roi[,c("Start","Stop","Chr","Length","ID", "Class")]
	
	tmp2 <- tmp[with(tmp, order(!Length, Start, Stop)), ]  # , -Length

	Samples <- 1:length(unique(tmp2$ID))
	names(Samples) <- unique(tmp2$ID)

	tmp2$Indx <- sapply(tmp2$ID, function(X){ Samples[X] })
	CytoBands$Indx <- rep(-1, nrow(CytoBands))
	roi$Indx <- rep(-2, nrow(roi))
	
	tmp2 <- rbind(tmp2, CytoBands, roi)
	
	#for(i in 1:22){ tmp2$Indx[tmp2$Chr %in% i & tmp2$Indx > 0] <- 1:length(tmp2$Indx[tmp2$Chr %in% i & tmp2$Indx > 0]) }

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

	tmp2ROI <- subset(tmp2, Class %in% "ROI")
	
	Colors <- brewer.pal("Set1", n=8) # + scale_color_brewer(palette="Set1")
	b <- ggplot(tmp2, aes(Start, Indx))
	b <- b + geom_segment(aes(x = Start, y = Indx, xend = Stop, yend = Indx, colour=as.factor(Class))) 
	b <- b + scale_colour_manual(values = c("ROI" = Colors[6],"q" = Colors[5],"p" = Colors[4], "1" = Colors[1], "3" = Colors[2], "4" = Colors[3], "0"=Colors[7], "2"=Colors[8]))
	b <- b + facet_wrap(~ Titles, scales = "free", ncol = NCOL) 
	b <- b + geom_vline(aes(xintercept = c(Start, Stop)), tmp2ROI, alpha=0.2) 
	
	ggsave(b, file=Name, width=width, height=height, dpi=300)
	return(tmp2)
}

		 
	

	
