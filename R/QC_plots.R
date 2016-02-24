##' QC_plots: QC plots. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title QC_plots
##' @param tmp: Unknown, default = Unknown. 
##' @return Unknown.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

QC_Plots <- function(tmp)
{
	library(gplots)
	Cor <- cor(tmp, method="s", use="p")
	tmp2 <- apply(Cor, 1, function(X){ names(X[abs(X) > 0.7 & abs(X) < 0.99]) })

	png("Heatmap400Samples.png", height=800, width=800)
	heatmap.2(Cor, trace="none", col=redblue, margins=c(11,11))
	dev.off()
	ColorsChr <- rainbow(22, alpha=0.7)

 
	names(ColorsChr) <- 1:22
	ColorsVector <- sapply(tmp$Chr, function(X){ ColorsChr[X] })
	tmp$ColorChr <- ColorsVector



	png("PointTestxLRR_SD.png", height=800, width=800, res=120)
	plot(tmp$Log.R.Ratio.SD, log10(tmp$PointTest)*-1, type="n")
	text(tmp$Log.R.Ratio.SD, log10(tmp$PointTest)*-1, labels=tmp$Chr, col=tmp$ColorChr, pch=1,cex=1)
	dev.off()
	cat("PointTestxLRR_SD\n")
	
	png("BBBBxLRR_Median.png", height=800, width=800, res=120)
	plot(tmp$BBBB, tmp$Log.R.Ratio.Median, col=tmp$ColorChr, type="n")
	text(tmp$BBBB, tmp$Log.R.Ratio.Median, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()
	cat("BBBBxLRR_Median\n")

	png("BBBBxBAF_Mean.png", height=800, width=800, res=120)
	plot(tmp$BBBB, tmp$B.Allele.Freq.Mean, col=tmp$ColorChr, type="n")
	text(tmp$BBBB, tmp$B.Allele.Freq.Mean, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()
	cat("BBBBxBAF_Mean.png\n")

	png("MeanGCxLRR_Mean.png", height=800, width=800, res=120)
	plot(tmp$MeanGC, tmp$Log.R.Ratio.Mean, col=tmp$ColorChr, type="n")
	text(tmp$MeanGC, tmp$Log.R.Ratio.Mean, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()
	cat("BBBBxBAF_Mean.png\n")

	png("MeanGCxBAF_Mean.png", height=800, width=800, res=120)
	plot(tmp$MeanGC, tmp$B.Allele.Freq.Mean, col=tmp$ColorChr, type="n")
	text(tmp$MeanGC, tmp$B.Allele.Freq.Mean, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("MeanGCxAmplitude.SD.png", height=800, width=800, res=120)
	plot(tmp$MeanGC, tmp$Amplitude.SD, col=tmp$ColorChr, type="n")
	text(tmp$MeanGC, tmp$Amplitude.SD, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("MeanGCxBBBB.png", height=800, width=800, res=120)
	plot(tmp$MeanGC, tmp$BBBB, col=tmp$ColorChr, type="n")
	text(tmp$MeanGC, tmp$BBBB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("MeanGCxAB.png", height=800, width=800, res=120)
	plot(tmp$MeanGC, tmp$AB, col=tmp$ColorChr, type="n")
	text(tmp$MeanGC, tmp$AB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("LRR_SDxBAF_Drift_SD.png", height=800, width=800, res=120)
	plot(tmp$Log.R.Ratio.SD, tmp$BAF_Drift_SD, col=tmp$ColorChr, type="n")
	text(tmp$Log.R.Ratio.SD, tmp$BAF_Drift_SD, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("LRR_SDxAAAB.png", height=800, width=800, res=120)
	plot(tmp$Log.R.Ratio.SD, tmp$AAAB, col=tmp$ColorChr, type="n")
	text(tmp$Log.R.Ratio.SD, tmp$AAAB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("BAF_MedianxAAAA.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Median, tmp$AAAA, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Median, tmp$AAAA, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("BAF_MedianxAB.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Median, tmp$AB, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Median, tmp$AB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("BAF_MedianxBBBB.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Median, tmp$BBBB, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Median, tmp$BBBB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	# Mean
	png("BAF_MeanxAAAA.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Mean, tmp$AAAA, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Mean, tmp$AAAA, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()
	
	png("BAF_MeanxAB.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Mean, tmp$AB, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Mean, tmp$AB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()

	png("BAF_MeanxBBBB.png", height=800, width=800, res=120)
	plot(tmp$B.Allele.Freq.Mean, tmp$BBBB, col=tmp$ColorChr, type="n")
	text(tmp$B.Allele.Freq.Mean, tmp$BBBB, col=tmp$ColorChr, labels=tmp$Chr)
	dev.off()
	
}
