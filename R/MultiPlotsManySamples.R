MultiPlotsManySamples <- function(df=iPsych_8_Samples_Good, Name="Test.png")
{
	library(ggplot2)
	library(RColorBrewer)
	df <- df[,c("Start", "Stop", "Chr", "ID", "CN")]
	Cyto <- CytoBands
	Cyto$CN <- Cyto$Class
	Cyto <- Cyto[,c("Start", "Stop", "Chr", "ID", "CN")]
	Cyto$Chr2 <- factor(Cyto$Chr, levels=1:22)
	df$Chr2 <- factor(df$Chr, levels=1:22)

	Colors <- brewer.pal("Set1", n=9) # + scale_color_brewer(palette="Set1")

	plots <- list()  # new empty list
	IDs <- sort(unique(df$ID))
	for (i in 1:length(unique(df$ID))) 
	{
		tmp <- subset(df, ID %in% IDs[i])
		p1 <- ggplot(Cyto) + geom_segment(aes(x = Start, y = 0, xend = Stop, yend = 0, colour = as.factor(CN))) + facet_wrap(~Chr2, nrow=1, scales = "free") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position="none") + ylim(0,1) + ylab(IDs[i]) + xlab("")  + scale_colour_manual(values = c("q" = Colors[5],"p" = Colors[4], "1" = Colors[1], "3" = Colors[2], "4" = Colors[3], "0"=Colors[7], "2"=Colors[9]))  #+ ggtitle(IDs[i])
		p1 <- p1 + geom_vline(aes(xintercept = c(Start), colour=factor(CN)), data=tmp)
	   	plots[[i]] <- p1  # add each plot into plot list
	}
	png(Name, height=1200, width=1800, res=100)
	multiplot(plotlist = plots, cols = 1)
	dev.off()
}
