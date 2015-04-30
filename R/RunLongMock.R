##' RunLongMock: Run a single mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RunSingleMock
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

RunLongMock <- function(Noise=0, CNVMean=0.4, Name="Test", Method="PennCNV", HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/" )
{
	LongRoi <- MakeLongMockSample(Size=500)
	Sample <- read.table("LongMockSample.tab", sep="\t", header=TRUE, stringsAsFactors=F)
	
	if(Method %in% "PennCNV")
	{
		PredictedCNV <- RunPennCNV(PathRawData=".", Pattern="LongMockSample.tab$", Cores=1, Skip=0, Normalization=FALSE, PFB="")
	}
	else if(Method %in% "iPsychCNV")
	{
		PredictedCNV <- iPsychCNV(PathRawData=".", MINNumSNPs=28, Cores=1, Pattern="LongMockSample.tab$", MinLength=10, Skip=0, LCR=FALSE, Path2PennCNV=Path2PennCNV, HMM=HMM )
	}
	else if(Method %in% "Gada")
	{
		PredictedCNV <- RunGada(Sample)
	}	

	Name <- paste(Name, "_Noise", Noise, "_CNVmean", CNVMean, ".png", sep="", collapse="")
	PlotLRRAndCNVs(PredictedCNV, Sample, CNVMean, Name=Name, Roi=LongRoi)
}
