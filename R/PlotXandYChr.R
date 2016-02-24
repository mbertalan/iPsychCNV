##' PlotXandYChr: Plot Log R Ratio (LRR) and copy number variations (CNVs) for the X and Y chromosome the whole mock sample.   
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title PlotXandYChr
##' @param Sample: Unknown, default = Unknown.
##' @param Name: Unknown, default = Unknown.
##' @return png plot
##' @author Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples unknown
##' 

PlotXandYChr <- function(Sample=MockData, Name="Test.png")
{
	library(iPsychCNV)
	library(ggplot2)
	library(ggbio)
	library(RColorBrewer)

	tmp <- subset(Sample, Chr %in% "X")
	tmp2 <- subset(Sample, Chr %in% "Y")

	tmp <- tmp[order(tmp$Position),]
	tmp2 <- tmp2[order(tmp2$Position),]

	# Fixing Log.R.Ratio for crazy values (deCODE with LRR == 999)
	tmp$Log.R.Ratio[abs(tmp$Log.R.Ratio) > 10] <- 0
	tmp2$Log.R.Ratio[abs(tmp2$Log.R.Ratio) > 10] <- 0

	Colors <- brewer.pal("Set1", n=9)

	Mean <- SlideWindowMean(tmp$Log.R.Ratio, 50)
	tmp$Mean <- Mean

	Mean <- SlideWindowMean(tmp2$Log.R.Ratio, 25)
	tmp2$Mean <- Mean

	## X ##
	# CHR
	chr <- unique(tmp$Chr)
	data(hg19IdeogramCyto,package="biovizBase")
	CHR <- paste("chr", chr, collapse="", sep="")
	p3 <- plotIdeogram(hg19IdeogramCyto,CHR,cytoband=TRUE,xlabel=TRUE, aspect.ratio = 1/85, alpha = 0, color="gray") # alpha = 0, color="black" 

	# Limits
	Min <- median(tmp$Log.R.Ratio) - (sd(tmp$Log.R.Ratio)*4)
	Max <- median(tmp$Log.R.Ratio) + (sd(tmp$Log.R.Ratio)*4)
	#Min <- median(tmp$Log.R.Ratio) - (quantile(tmp$Log.R.Ratio, probs=2/100)*2)
	#Max <- median(tmp$Log.R.Ratio) + (quantile(tmp$Log.R.Ratio, probs=98/100)*2)

	# LRR
	p1 <- ggplot(tmp, aes(Position, Log.R.Ratio))
	p1 <- p1 + geom_point(alpha=0.05, size=1, colour=Colors[1]) + ylim(Min, Max)	
	p1 <- p1 + geom_line(data=tmp, aes(x=Position, y = Mean), size = 0.5, alpha=0.5, colour="black") 

	# BAF
	p2 <- ggplot(tmp, aes(Position, y = B.Allele.Freq)) 
	p2 <- p2 + geom_point(aes(col="B.Allele.Freq"), alpha=0.3, size=1, colour=Colors[2])  
	p2 <- p2 + scale_colour_manual(values = c("B.Allele.Freq"=Colors[2], "Mean"="black", Log.R.Ratio=Colors[1]))		
	

	## Y ##
	# CHR
	chr <- unique(tmp2$Chr)
	data(hg19IdeogramCyto,package="biovizBase")
	CHR <- paste("chr", chr, collapse="", sep="")
	p4 <- plotIdeogram(hg19IdeogramCyto,CHR,cytoband=TRUE,xlabel=TRUE, aspect.ratio = 1/85, alpha = 0, color="gray") # alpha = 0, color="black" 

	# Limits
	Min <- median(tmp2$Log.R.Ratio) - (sd(tmp2$Log.R.Ratio)*4)
	Max <- median(tmp2$Log.R.Ratio) + (sd(tmp2$Log.R.Ratio)*4)
	#Min <- median(tmp$Log.R.Ratio) - (quantile(tmp$Log.R.Ratio, probs=2/100)*2)
	#Max <- median(tmp$Log.R.Ratio) + (quantile(tmp$Log.R.Ratio, probs=98/100)*2)

	# LRR
	p5 <- ggplot(tmp2, aes(Position, Log.R.Ratio))
	p5 <- p5 + geom_point(alpha=0.05, size=1, colour=Colors[1]) + ylim(Min, Max)		
	p5 <- p5 + geom_line(data=tmp2, aes(x=Position, y = Mean), size = 0.5, alpha=0.5, colour="black") 
	
	# BAF
	p6 <- ggplot(tmp2, aes(Position, y = B.Allele.Freq)) 
	p6 <- p6 + geom_point(aes(col="B.Allele.Freq"), alpha=0.5, size=1, colour=Colors[2])  
	p6 <- p6 + scale_colour_manual(values = c("B.Allele.Freq"=Colors[2], "Mean"="black", Log.R.Ratio=Colors[1]))
	
	# Sample Name
	ID <- unique(Sample$Sample.ID)[1]
	Title <- paste("Sample: ", ID, sep="", collapse="")

	Plot <- tracks(p3, p1, p2, p4, p5, p6, title=Title)

	Plot@fixed <- rep(TRUE, 6)
	Plot@hasAxis <- c(TRUE, TRUE,TRUE,TRUE,TRUE, FALSE)
	
	ggsave(Name, plot=Plot, height=8, width=16, dpi=100) 
}
