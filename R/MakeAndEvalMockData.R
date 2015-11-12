##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

MakeAndEvalMockData <- function(N=100, Wave1=FALSE, Type="PKU", Cores=20)
{
	library(iPsychCNV)
	# Creating Mock Data
	MockFileName <- paste(Type, "_Mock_CNVs", ".RData", sep="", collapse="")
	Mock_CNVs <- MockData(N, Wave1, Type, Cores)
	save(Mock_CNVs, file=MockFileName)

	# Creating Roi from Mock Data
	Mock_CNV.Roi <- subset(Mock_CNVs, ID %in% "MockSample_1.tab" & CN != 2)
	Mock_CNV.Roi$Class <- rep("ROI", nrow(Mock_CNV.Roi))
	save(Mock_CNV.Roi, file="Mock_CNV.Roi.RData")

	# QC
	QC_Mock <- QualityControl(PathRawData=".", Cores=Cores, Pattern="MockSample.*\\.tab$") 
	MockFileName <- paste("QC_", Type, "_Mock", ".RData", sep="", collapse="")
	TypeName <- paste(Type, "_Mock", sep="", collapse="")
	QC_Mock$Source <- rep(TypeName, nrow(QC_Mock))
	QC_Mock$Source2 <- rep(Type, nrow(QC_Mock))
	save(QC_Mock, file=MockFileName)
	

	# Running penncnv
	PennCNV_Mock_CNVs <- RunPennCNV(PathRawData=".", Pattern="MockSample.*\\.tab$", Cores=Cores, Skip=0, Normalization=FALSE, PFB="")
	PennFileName <- paste("PennCNV_", Type, "_Mock_CNVs", ".RData", sep="", collapse="")
	save(PennCNV_Mock_CNVs, file=PennFileName)
	
	# Plot all CNVs
	PlotFileName <- paste("PennCNV_", Type, "_Mock_CNVs", ".png", sep="", collapse="")
	PlotAllCNVs(PennCNV_Mock_CNVs, Name=PlotFileName, Roi=Mock_CNV.Roi)

	# Evaluating Program result
	PennCNV_Mock_CNVs.Eval <- EvaluateMockResults(Mock_CNVs, PennCNV_Mock_CNVs)
	
	save(PennCNV_Mock_CNVs.Eval, file="PennCNV_Mock_CNVs.Eval")

	# ROC
	RocFileName <- paste("ROC_PennCNV_", Type, "_Mock_CNVs", ".png", sep="", collapse="")
	png(RocFileName, height=800, width=800)
	plot.roc(as.factor(PennCNV_Mock_CNVs.Eval$CNV.Present), PennCNV_Mock_CNVs.Eval$Overlap.SNP, print.auc=TRUE, percentage=TRUE)
	dev.off()
}
