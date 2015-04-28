##' RunSingleMock: Run a single mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RunSingleMock
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

RunSingleMock <- function(Noise=0, CNVMean=0.4, Name="Test", Method="PennCNV")
{
	tmp <- MakeOneMockSample(Noise=Noise, CNVMean=CNVMean)
	if(Method %in% "PennCNV")
	{
		PredictedCNV <- RunPennCNV(PathRawData=".", Pattern="MockSample_1.tab$", Cores=1, Skip=0, Normalization=FALSE, PFB="")
	}
	else if(Method %in% "iPsychCNV")
	{
		PredictedCNV <- iPsychCNV(PathRawData=".", MINNumSNPs=28, Cores=1, Pattern="MockSample_1.tab$", MinLength=10, Skip=0, LCR=FALSE)
	}
	else if(Method %in% "Gada")
	{
		PredictedCNV <- RunGada(tmp)
	}	

	Name <- paste(Name, "_Noise", Noise, "_CNVmean", CNVMean, ".png", sep="", collapse="")
	PlotLRRAndCNVs(PredictedCNV, tmp, CNVMean, Name=Name, Roi=RoiSingleMock)
}
