##' MultiPlotsManySamples: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MultiPlotsManySamples
##' @param df: Data frame with predicted CNVs for each sample, default = Unknown.
##' @param Name: Unknown, default =Unknown.
##' @param Height: Unknown, default = 1800.
##' @param Width: Unknown, default = 2000.
##' @param Res: Resolution, defualt = 120.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

MultiPlotsManySamples <- function(df=iPsych_8_Samples_Good, Name="Test.png", height=1800, width=2000, res=120)
{
	library(ggplot2)
	library(RColorBrewer)
	df <- df[,c("Start", "Stop", "Chr", "ID", "CN", "Source")]
	df$Chr2 <- factor(df$Chr, levels=1:22)
	
	# Cytobands
	Cyto <- CytoBands
	Cyto$CN <- Cyto$Class
	
	#Including source on Cyto
	CytoAll <- sapply(unique(df$Source), function(X){ Cyto$Source <- X; return(Cyto) })
	Cyto <- MatrixOrList2df(CytoAll)

	Cyto <- Cyto[,c("Start", "Stop", "Chr", "ID", "CN", "Source")]
	Cyto$Chr2 <- factor(Cyto$Chr, levels=1:22)
	Cyto <- subset(Cyto, Source %in% unique(df$Source))
	
	Colors <- brewer.pal("Set1", n=9) # + scale_color_brewer(palette="Set1")

	plots <- list()  # new empty list
	IDs <- sort(unique(df$ID))
	for (i in 1:length(unique(df$ID))) 
	{
		tmp <- subset(df, ID %in% IDs[i])
		#p1 <- ggplot() + geom_segment(data=Cyto, aes(x = Start, y = 0, xend = Stop, yend = 0, colour = as.factor(CN))) + facet_wrap(Chr2~Source, nrow=1, scales = "free") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position="none") + ylim(0,1) + ylab(IDs[i]) + xlab("")  + scale_colour_manual(values = c("q" = Colors[5],"p" = Colors[4], "1" = Colors[1], "3" = Colors[2], "4" = Colors[3], "0"=Colors[7], "2"=Colors[9]))  #+ ggtitle(IDs[i])
		#p1 <- p1 + geom_vline(aes(xintercept = c(Start), colour=factor(CN)), data=tmp)
		p1 <- ggplot() + geom_segment(data=Cyto, aes(x = Start, y = 0, xend = Stop, yend = 0, colour = as.factor(CN))) + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position="none") + ylim(0,0.2) + ylab(IDs[i]) + xlab("")  + scale_colour_manual(values = c("q" = Colors[5],"p" = Colors[4], "1" = Colors[1], "3" = Colors[2], "4" = Colors[3], "0"=Colors[7], "2"=Colors[9]))  #+ ggtitle(IDs[i])
		p1 <- p1 + geom_vline(aes(xintercept = c(Start), colour=factor(CN)), data=tmp) + facet_grid(Source~Chr2, scales="free")
		
	   	plots[[i]] <- p1  # add each plot into plot list
	}
	png(Name, height=height, width=width, res=res)
	Multiplot(plotlist = plots, cols = 1)
	dev.off()
}
