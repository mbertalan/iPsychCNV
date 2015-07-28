##' ReScanCNVs: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ReScanCNVs
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

ReScanCNVs <- function(CNVs=CNVs, PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", MINNumSNPs=20, Cores=1, hg="hg18", NumFiles="All", Pattern="22q11_*", MinLength=10, SelectedFiles=NA, Skip=10, LCR=FALSE, PFB=NULL, chr=NA, penalty=60, Quantile=FALSE, QSpline=FALSE, sd=0.18, recursive=FALSE, CPTmethod="meanvar", CNVSignal=0.1, penvalue=10, OutputPath=NA) # Files2 OutputPath
{	
	if(file.exists("Progress.txt")){ file.remove("Progress.txt") }
	suppressPackageStartupMessages(library(parallel))
	ptm <- proc.time()
	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
	
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
		CNVs <- GetIndxPositionFromChips(CNVs, Sample)
		
		if(nrow(CNVs) > 0)
		{
			CNVs <- subset(CNVs, NumSNPs > MINNumSNPs)
			df <- FilterCNVs.V4(CNVs = CNVs, MinNumSNPs=MINNumSNPs, Sample=Sample, ID) # PathRawData = PathRawData,
			df$Source <- rep("iPsychCNV", nrow(df))
			df$CN <- df$Class
			df$CN[df$CN %in% "Del"] <- "1"
			df$CN[df$CN %in% "Normal"] <- "2"
			df$CN[df$CN %in% "Dup"] <- "3"
			df$CN[df$CN %in% "DoubleDel"] <- "0"
			df$CN[df$CN %in% "DoubleDup"] <- "4"
			df$CN <- as.numeric(df$CN)
			df$SampleID <- ID
			return(df)
		}	
		Res.tmp <- proc.time() - ptm.tmp
	})
	cat("Done all !\n")
	df <- MatrixOrList2df(tmp)
	df <- subset(df, CN != 2) # removing non-CNV results to save memory
	df <- df[, !colnames(df) %in% "Class"]
	
	if(!is.na(OutputPath))
	{
		# Print results or return as an object ?
		OutputFile <- paste(OutputPath, ID, ".CNVs", sep="", collapse="")
		write.table(df, file=OutputFile, sep="\t", quote=FALSE, row.names=FALSE)
	}
		
	TimeRes <- proc.time() - ptm
	cat("Total time: ", TimeRes["elapsed"], "\n")
	return(df)
}
