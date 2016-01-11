##' ReScanCNVs: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ReScanCNVs
##' @param CNVs data frame with hotspots to re-scan. Minimum information necessary is Chr, Start and Stop position. 
##' @param PathRawData The path to the raw data files contining LRR and BAF values.
##' @param MINNumSNPs Minimum number of SNPs per CNV. Default 20.
##' @param Cores Number of cores used. Default 1.
##' @param hg Human genome version. Default hg19.
##' @param NumFiles Number of files to be analyzed from PathRawData.
##' @param Pattern File pattern in the PathRawData. Example: "*.txt".
##' @param MinLength Minimum CNV length.
##' @param SelectedFiles List of file names that should be analyzed from PathRawData. 
##' @param Skip integer: the number of lines of the data file to skip before beginning to read data.
##' @param LCR list: low complex region, list of SNPs that should be removed.
##' @param PFB vector: Population frequency 0 to 1 for each SNP in the array. 
##' @param chr character: select a specific chromosome to be analyzed. 
##' @param penalty the coefficient of the penalty for degrees of freedom in the GCV criterion. From smooth.spline {stats}.
##' @param Quantile logical, if quantile normalization should be applied or not. Default FALSE.
##' @param QSpline logical, if a cubic smoothing spline should be used to normalize the data. Default FALSE.
##' @param sd numeric Log R ratio standard deviation for the quantile nomarlization. Default 0.18.
##' @param recursive logical, Should the listing recurse into directories? From list.files {base}.
##' @param CPTmethod character, method to find change points from changepoint package by Rebecca Killick. Default "meanvar", or "mean".
##' @param CNVSignal numeric, minumum CNV signal to be consider a CNV in absolute value. Default 0.1, any CNV with mean Log R ration in the CNV region with abs(X) < 0.1 is ignored. 
##' @param penvalue Same as pen.value from function cpt.mean at changepoint R package by Rebecca Killick. Default 10. "The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  A vector of length 2 (min,max) if using the CROPS penalty.  The value of the penalty when using the Manual penalty option - this can be a numeric value or text           giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters".
##' @param OutputPath character, path for output.
##' @param OutputFileName character, Output file name. 
##' @param OnlyCNVs logical, if TRUE only CNVs with copy number state 0,1,3,4 will be returned. If FALSE will return also changepoint regions with CN = 2.
##' @param IndxPos, if index position for each hotspot is not included it will calculate. However, this step is time consuming.
##' @param ResPerSample, If TRUE saves the results from each sample.  
##' @param Files, a vector with all samples name and path. If too many samples, list all files using recursive=TRUE can take long time. 
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' MockDataCNVs <- MockData(N=100, Type="PKU", Cores=20)
##' iPsych.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=20, Cores=20, Pattern="^MockSample", MinLength=10, Skip=0)
##' iPsych.Pred.hotspots <- HotspotsCNV(iPsych.Pred, Freq=2, OverlapCutoff=0.9, Cores=1)
##' iPsych.Pred.Rescan <- ReScanCNVs(iPsych.Pred.hotspots, PathRawData="./Data", hg="hg19", Pattern="*", Cores=20, Skip=0)
##' iPsych.Pred.Rescan$ID <- iPsych.Pred.Rescan$SampleID
##' PlotAllCNVs(iPsych.Pred.Rescan, Name="iPsych.Pred.Rescan.png", hg="hg19", Roi=MockDataCNVs.roi)

ReScanCNVs <- function(CNVs=CNVs, PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", MINNumSNPs=5, Cores=1, hg="hg19", NumFiles="All", Pattern="*", MinLength=10, SelectedFiles=NA, Skip=10, LCR=FALSE, PFB=NULL, chr=NA, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18, recursive=FALSE, CPTmethod="meanvar", CNVSignal=0.1, penvalue=10, OutputPath=NA, IndxPos=FALSE, ResPerSample=FALSE, Files=NA) # Files2 OutputPath
{	
	if(file.exists("Progress.txt")){ file.remove("Progress.txt") }
	suppressPackageStartupMessages(library(parallel))
	ptm <- proc.time()
	
	if(is.na(Files))
	{
		Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
	}
	
	if(length(SelectedFiles) > 1 & !is.na(SelectedFiles[1]))
	{
		Files <- sapply(SelectedFiles, function(X){ Files[grep(X, Files)] })
	}

	if(NumFiles %in% "All"){ NumFiles <- length(Files) }

	cat("Running ", NumFiles, "files\n")
	tmp <- mclapply(Files[1:NumFiles], mc.cores=Cores, mc.preschedule = FALSE, function(RawFile) 
	{
		write(RawFile,file="Progress.txt",append=TRUE)
		Count <- length(readLines("Progress.txt"))	
		Percent <- round((Count/NumFiles)*100)
		Percent <- paste(Percent, "%", sep="", collapse="")
		cat("Running:\t", RawFile, "\t\t", Percent, "\n")
		ID <- tail(unlist(strsplit(RawFile, "/")),n=1)
	
		# Read sample file		
		ptm.tmp <- proc.time()
		Sample <- ReadSample(RawFile, skip=Skip, LCR=LCR, PFB=PFB, chr=chr)

		# The CNVs are given by hotspots, Position might change if multiple chips are used.
		if(IndxPos)
		{
			CNVs <- GetIndxPositionFromChips(CNVs, Sample)
		}
		
		if(nrow(CNVs) > 0)
		{
			CNVs <- subset(CNVs, NumSNPs > MINNumSNPs)
			df <- FilterCNVs.V4(CNVs = CNVs, MinNumSNPs=MINNumSNPs, Sample=Sample, ID) # PathRawData = PathRawData,

			# removing non-CNV results to save memory
			df <- subset(df, CN != 2) 
			
			# Save memory if too many CNVRs to ReScan. Print results or return as an object ?
			if(ResPerSample)
			{
				OutputFile <- paste(OutputPath, ID, ".CNVs", sep="", collapse="")
				write.table(df, file=OutputFile, sep="\t", quote=FALSE, row.names=FALSE)
			}
			else
			{
				return(df)
			}
		}	
		Res.tmp <- proc.time() - ptm.tmp
	})
	cat("Done all !\n")
	df <- MatrixOrList2df(tmp)
	
	if(!is.na(OutputPath))
	{
		OutputFile <- paste(OutputPath, "All_CNVs.tab", sep="", collapse="")
		write.table(df, file=OutputFile, sep="\t", quote=FALSE, row.names=FALSE)
	}
		
	TimeRes <- proc.time() - ptm
	cat("Total time: ", TimeRes["elapsed"], "\n")
	return(df)
}
