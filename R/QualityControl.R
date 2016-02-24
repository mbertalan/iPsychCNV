##' QualityControl: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title QualityControl
##' @param PathRawData: The path to the raw data files containing the Log R Ratio (LRR) and B Allele Frequency (BAF) vaules. 
##' @param Cores: Number of cores used, default = 10.
##' @param Pattern: File pattern in the raw data, default = "*_*assed_*".
##' @param NumFiles: Number of files to be analyzed from raw data, default = All.
##' @param Skip: Integer, the number of lines of the data file to be skipped before beginning to read the data, default = 0.
##' @param Normalization: Unknown, default = FALSE. 
##' @return Data frame with QC variables.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

QualityControl <- function(PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", Cores=10, Pattern="*_*assed_*", NumFiles="All", Skip=0, Normalization=FALSE) # Files2 OutputPath
{	
	library(parallel)

	Files <- list.files(PathRawData, pattern=Pattern, recursive=T, full.names=T) #  *_*assed_*

	if(NumFiles %in% "All"){ NumFiles <- length(Files) }


	tmp <- mclapply(Files[1:NumFiles], mc.cores=Cores, mc.preschedule = FALSE, function(X) 
	{
		
		ID <- tail(unlist(strsplit(X, "/")),n=1)
		cat(ID,"\n")
	
		CNV <- ReadSample(X, skip=Skip)
		
		if(Normalization)
		{
			CNV <- NormalizeData(CNV, ExpectedMean=0, DF=NA, NormQspline=FALSE)
		}

		Res <- PeaksAndPits(CNV, ID)
		return(Res)
	})
	tmp2 <- MatrixOrList2df(tmp)
	tmp3 <- tmp2[,!colnames(tmp2) %in% ".id"]
	indx <- which(is.na(as.numeric(tmp3$Log.R.Ratio.Mean)))
	if(length(indx) > 0)
	{
		tmp3 <- tmp3[indx*-1,]
		tmp <- sapply(tmp3[,!colnames(tmp3) %in% "ID"], as.numeric)
		tmp <- as.data.frame(tmp)
		tmp$ID <- tmp3$ID
		tmp3 <- tmp
	}
	return(tmp3)
}
