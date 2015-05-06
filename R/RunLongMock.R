##' RunLongMock: Run a single mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RunSingleMock
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

RunLongMock <- function(Name="Test", Method="PennCNV", CNVDistance=1000, Type=c(0,1,2,3,4), Mean=c(-0.3, -0.6, 0.3, 0.6), Size=c(300, 600), HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/" )
{
	LongRoi <- MakeLongMockSample(CNVDistance, Type, Mean, Size)
	Sample <- read.table("LongMockSample.tab", sep="\t", header=TRUE, stringsAsFactors=F)
	
	if(Method %in% "PennCNV")
	{
		PredictedCNV <- RunPennCNV(PathRawData=".", Pattern="^LongMockSample.tab$", Cores=1, Skip=0, Normalization=FALSE, PFB=NA, HMM=HMM, Path2PennCNV=Path2PennCNV )
	}
	else if(Method %in% "iPsychCNV")
	{
		PredictedCNV <- iPsychCNV(PathRawData=".", MINNumSNPs=28, Cores=1, Pattern="^LongMockSample.tab$", MinLength=10, Skip=0, LCR=FALSE, Quantile=FALSE, NormMean=FALSE)
	}
	else if(Method %in% "Gada")
	{
		PredictedCNV <- RunGada(Sample)
	}	

	Name <- paste(Name,"_",Method, "_LongMockResult.png", sep="", collapse="")
	CNVMean=0.3
	PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
}
