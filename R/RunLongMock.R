##' RunLongMock: Run a single mock sample.   
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RunSingleMock
##' @param Name: Unknown.
##' @param Method: The Copy Number Variation (CNV) calling algorithm, default = PennCNV.
##' @param CNVDistance: Distance between CNVs, default = 1000.
##' @param Type: Unknown, default = Unknown.
##' @param Mean: Unknown, default = Unknown.
##' @param Size: Unknown, default = Unknown.  
##' @param HMM: Unknown, default = Unknown. 
##' @param Path2PennCNV: The path to the pennCNV algorithm, Unknown. 
##' @return Classification for Log R Ratio (LRR).
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

RunLongMock <- function(Name="Test", Method="PennCNV", CNVDistance=1000, Type=c(0,1,2,3,4), Mean=c(-0.3, -0.6, 0.3, 0.6), Size=c(300, 600), HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/" )
{
	LongRoi <- MakeLongMockSample(CNVDistance, Type, Mean, Size)
	Sample <- read.table("LongMockSample.tab", sep="\t", header=TRUE, stringsAsFactors=F)
	Name <- paste(Name,"_",Method, "_LongMockResult.png", sep="", collapse="")
	CNVMean=0.3
	
	if(Method %in% "PennCNV")
	{
		PredictedCNV <- RunPennCNV(PathRawData=".", Pattern="^LongMockSample.tab$", Cores=1, Skip=0, Normalization=FALSE, PFB=NA, HMM=HMM, Path2PennCNV=Path2PennCNV )
		PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
	}
	else if(Method %in% "iPsychCNV")
	{
		PredictedCNV <- iPsychCNV(PathRawData=".", MINNumSNPs=28, Cores=1, Pattern="^LongMockSample.tab$", MinLength=10, Skip=0, LCR=FALSE, Quantile=FALSE)
		PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
	}
	else if(Method %in% "Gada")
	{
		PredictedCNV <- RunGada(Sample)
		PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
	}	
	else if(Method %in% "All")
	{
		PennCNV.Pred <- RunPennCNV(PathRawData=".", Pattern="^LongMockSample.tab$", Cores=1, Skip=0, Normalization=FALSE, PFB=NA, HMM=HMM, Path2PennCNV=Path2PennCNV )
		PlotLRRAndCNVs(PennCNV.Pred, Sample, CNVMean, Name="PennCNV.Pred.png", Roi=LongRoi)
		PennCNV.Eval <- EvaluateMockResults(LongRoi, PennCNV.Pred)
		
		iPsychCNV.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=28, Cores=1, Pattern="^LongMockSample.tab$", MinLength=10, Skip=0, LCR=FALSE, Quantile=FALSE)
		PlotLRRAndCNVs(iPsychCNV.Pred, Sample, CNVMean, Name="iPsychCNV.Pred.png", Roi=LongRoi)
		iPsychCNV.Eval <- EvaluateMockResults(LongRoi, iPsychCNV.Pred)
		
		Gada.Pred <- RunGada(Sample)
		Gada.Pred$ID <- "LongMockSample.tab"
		PlotLRRAndCNVs(Gada.Pred, Sample, CNVMean, Name="Gada.Pred.png", Roi=LongRoi)
		Gada.Eval <- EvaluateMockResults(LongRoi, Gada.Pred)
		
		ColNames <- intersect(intersect(colnames(iPsychCNV.Pred), colnames(PennCNV.Pred)), colnames(Gada.Pred))
		PredictedCNV <- rbind(PennCNV.Pred[,ColNames], iPsychCNV.Pred[,ColNames], Gada.Pred[,ColNames])
		#save(All.Pred, file="All.Pred.RData")
		png("ROC_All.png")
		plot.roc(PennCNV.Eval$CNV.Predicted, PennCNV.Eval$CNV.Present, percent=TRUE, col="#377eb8")
		plot.roc(iPsychCNV.Eval$CNV.Predicted, iPsychCNV.Eval$CNV.Present, percent=TRUE, col="#e41a1c", add=TRUE)
		plot.roc(Gada.Eval$CNV.Predicted, Gada.Eval$CNV.Present, percent=TRUE, col="#4daf4a", add=TRUE)
		legend("bottomright", legend=c("iPsychCNV", "PennCNV", "Gada"), col=c("#e41a1c","#377eb8", "#4daf4a"), lwd=2)
		dev.off()
	}
		

	#PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
	#CNVs.Eval <­ EvaluateMockResults(LongRoi, PredictedCNV)
	
	return(PredictedCNV)
}
