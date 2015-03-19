##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

iPsychCNV <- function(PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", MINNumSNPs=100, Cores=10, NumFiles="All", Pattern="22q11_*", MinLength=100000, SelectedFiles=NA, NormQspline=FALSE, DFspline=NA, Skip=10, LCR=TRUE) # Files2 OutputPath
{	
	if(file.exists("Progress.txt")){ file.remove("Progress.txt") }

	suppressPackageStartupMessages(library(parallel))
	
	ptm <- proc.time()

	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=TRUE)
	
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
	
		# Read CNV file		
		CNV <- ReadSample(RawFile, skip=Skip, LCR=TRUE)
		
		# Normalize data
		CNV <- NormalizeData(CNV, ExpectedMean=0, DF=DFspline, NormQspline)
		
		### FIND CNVs ###
		CNVs <- FindCNV.V4(ID, MINNumSNPs, CNV)
		# CNVs <- subset(CNVs, NumSNPs < 2000) # Also Max number of SNPs
		CNVs <- subset(CNVs, Length > MinLength)
		if(nrow(CNVs) > 0)
		{
			CNVsRes <- FilterCNVs.V4(CNVs = CNVs, MinNumSNPs=MINNumSNPs, CNV, ID) # PathRawData = PathRawData,
			return(CNVsRes)
		}
	})
	cat("Done all !\n")
	
	
	df <- MatrixOrList2df(tmp)

	### Re-Searching CNVs ###
	# df2 <- ReSearching(Files[1:NumFiles], PathRawData, Cores, df, MINNumSNPs)
	# df2$Steps <- rep("Second", nrow(df2))
	# save(df2, file="Df2.RData")

	df$Source <- rep("iPsychCNV", nrow(df))
	df$CN <- df$Class
	df$CN[df$CN %in% "Del"] <- "1"
	df$CN[df$CN %in% "Normal"] <- "2"
	df$CN[df$CN %in% "Dup"] <- "3"
	df$CN[df$CN %in% "DoubleDel"] <- "0"
	df$CN[df$CN %in% "DoubleDup"] <- "4"
		

	TimeRes <- proc.time() - ptm
	cat("Total time: ", TimeRes["elapsed"], "\n")
	return(df)
}
