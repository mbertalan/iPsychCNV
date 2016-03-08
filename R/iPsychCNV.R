##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays.
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @param PathRawData: The path to the raw data files contining Log R Ratio (LRR) and B Allele Frequency (BAF) values.
##' @param MINNumSNPs: Minimum number of SNPs per CNV, default = 20.
##' @param Cores: Number of cores used, default =  1.
##' @param Hg: Human genome version, default = hg19.
##' @param NumFiles: Number of files to be analyzed from PathRawData.
##' @param Pattern: File pattern in the PathRawData. Example: "*.txt".
##' @param MinLength: Minimum CNV length, default = Unknown.
##' @param SelectedFiles: List of file names that should be analyzed from PathRawData.
##' @param Skip: Integer, the number of lines of the data file to skip before beginning to read data.
##' @param LCR list: Low copy repeat region, list of SNPs that should be removed.
##' @param PFB vector: Population frequency 0 to 1 for each SNP in the array.
##' @param Chr: Character, select a specific chromosome to be analyzed.
##' @param Penalty: The coefficient of the penalty for degrees of freedom in the GCV criterion. From smooth.spline {stats}.
##' @param Quantile: Logical, if quantile normalization should be applied or not, default = FALSE.
##' @param QSpline: Logical, if a cubic smoothing spline should be used to normalize the data, default = FALSE.
##' @param Sd: Numeric, LRR standard deviation for the quantile nomarlization, default = 0.18.
##' @param Recursive: Logical, should the listing recurse into directories? From list.files {base}.
##' @param CPTmethod: Character, method to find change points from changepoint package by Rebecca Killick. Default "meanvar", or "mean".
##' @param CNVSignal: Numeric, minumum CNV signal to be consider a CNV in absolute value, default = 0.1, any CNV with mean LRR in the CNV region with abs(X) < 0.1 is ignored.
##' @param Penvalue: Same as pen.value from function cpt.mean at changepoint R package by Rebecca Killick, default = 10. "The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  A vector of length 2 (min,max) if using the CROPS penalty.  The value of the penalty when using the Manual penalty option - this can be a numeric value or text           giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters".
##' @param OutputPath: Character, path for output.
##' @param OutputFileName: Character, output file name.
##' @param OnlyCNVs: Logical, if TRUE only CNVs with copy number state 0,1,3,4 will be returned. If FALSE will return also changepoint regions with CN = 2.
##' @param SNPList: Getting Chr. and Position from another source than the RawFile - input should be the full path of the SNPList with columns: Name, Chr, amd Position. Any positions from the RawFile will be erased. A PFB-column is also allowed but will be overwritten by the PFB-parameter or exchanged with 0.5
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' mockCNV <- MockData(N=5, Type="Blood", Cores=1)
##' cnvs <- iPsychCNV(PathRawData=".", Cores=1, Pattern="^MockSample*", Skip=0)

iPsychCNV <- function(PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", MINNumSNPs=20, Cores=1, hg="hg19", NumFiles="All", Pattern="22q11_*", MinLength=10, SelectedFiles=NA, Skip=10, LCR=FALSE, PFB=NULL, chr=NA, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18, recursive=FALSE, CPTmethod="meanvar", CNVSignal=0.1, penvalue=20, OutputPath=NA, OutputFileName="Test", OnlyCNVs=FALSE, SNPList=NULL) # Files2 OutputPath
{
	if(file.exists("Progress.txt")){ file.remove("Progress.txt") }

	suppressPackageStartupMessages(library(parallel))

	ptm <- proc.time()

	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)

	if(length(SelectedFiles) > 1 & !is.na(SelectedFiles[1]))
	{
		Files <- sapply(SelectedFiles, function(X){ Files[grep(X, Files)] })
	}


	if(NumFiles %in% "All")
	{
		NumFiles <- length(Files)
	}

	cat("Running ", NumFiles, "files\n")
	tmp <- mclapply(Files[1:NumFiles], mc.cores=Cores, mc.preschedule = FALSE, function(X)
	{
		RawFile <- X
		write(X,file="Progress.txt",append=TRUE)
		Count <- length(readLines("Progress.txt"))
		Percent <- round((Count/NumFiles)*100)
		Percent <- paste(Percent, "%", sep="", collapse="")
		cat("Running:\t", X, "\t\t", Percent, "\n")
		ID <- tail(unlist(strsplit(X, "/")),n=1)

		# Read sample file
		ptm.tmp <- proc.time()
		Sample <- ReadSample(RawFile, skip=Skip, LCR=LCR, PFB=PFB, chr=chr, SNPList=SNPList)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Read Samples time: ", Res.tmp["elapsed"], "\n")

		#cat(nrow(Sample), "\n")

		# Normalize data
		ptm.tmp <- proc.time()
		Sample <- NormalizeData(Sample, ExpectedMean=0, penalty=penalty, Quantile=Quantile, QSpline=QSpline, sd=sd)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Normalization time: ", Res.tmp["elapsed"], "\n")

		#cat("Sample:", nrow(Sample), "\n")

		### FIND CNVs ###
		ptm.tmp <- proc.time()
		CNVs <- FindCNV.V4(ID=ID, Sample=Sample, CPTmethod=CPTmethod, CNVSignal=CNVSignal, penvalue=penvalue)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Find CNVs time: ", Res.tmp["elapsed"], "\n")

		#cat("CNVs:", nrow(CNVs), "\n")

		# Remove centromere
		CNVs <- RemoveCentromere(df=CNVs, HG=hg)

		#cat("CNVs:", nrow(CNVs), "\n")

		if(nrow(CNVs) > 0)
		{
			CNVs <- subset(CNVs, NumSNPs > MINNumSNPs)
			ptm.tmp <- proc.time()
			df <- FilterCNVs.V4(CNVs = CNVs, MinNumSNPs=MINNumSNPs, Sample=Sample, ID) # PathRawData = PathRawData,
			Res.tmp <- proc.time() - ptm.tmp

			df$Source <- rep("iPsychCNV", nrow(df))

			if(OnlyCNVs)
			{
				df <- subset(df, CN != 2) # removing non-CNV results to save memory
			}
			df <- df[, !colnames(df) %in% "Class"]

			# Check if every sample results should be printed or not.
			if(!is.na(OutputPath))
			{
				OutputFile <- paste(OutputPath, ID, ".CNVs", sep="", collapse="")
				write.table(df, file=OutputFile, sep="\t", quote=FALSE, row.names=FALSE)
			}
			return(df)
		}
	})
	cat("Done all !\n")
	if(length(tmp) == 0)
	{
		cat("Sorry no CNV found.\n")
		TimeRes <- proc.time() - ptm
		cat("Total time: ", TimeRes["elapsed"], "\n")
	}
	else
	{
		df <- MatrixOrList2df(tmp)

		if(!is.na(OutputPath))
		{
			write.table(df, file=OutputFile, sep="\t", quote=FALSE, row.names=FALSE)
			OutputFile <- paste(OutputPath,	OutputFileName, ".CNVs", sep="", collapse="")
		}
		TimeRes <- proc.time() - ptm
		cat("Total time: ", TimeRes["elapsed"], "\n")
		return(df)
	}
}
