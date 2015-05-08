##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export

iPsychCNV <- function(PathRawData = "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/", MINNumSNPs=100, Cores=10, NumFiles="All", Pattern="22q11_*", MinLength=100000, SelectedFiles=NA, Skip=10, LCR=TRUE, PFB=NULL, chr=NA, penalty=20, Quantile=TRUE, QSpline=TRUE, sd=0.2, recursive=FALSE) # Files2 OutputPath
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
		CNV <- ReadSample(RawFile, skip=Skip, LCR=LCR, PFB=PFB, chr=chr)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Read Samples time: ", Res.tmp["elapsed"], "\n")
		
		# Normalize data
		ptm.tmp <- proc.time()
		CNV <- NormalizeData(CNV, ExpectedMean=0, penalty=penalty, Quantile=Quantile, QSpline=QSpline, sd=sd)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Normalization time: ", Res.tmp["elapsed"], "\n")
		
		### FIND CNVs ###
		ptm.tmp <- proc.time()
		CNVs <- FindCNV.V4(ID, MINNumSNPs, CNV)
		Res.tmp <- proc.time() - ptm.tmp
		#cat("Find CNVs time: ", Res.tmp["elapsed"], "\n")
	
		CNVs <- subset(CNVs, Length > MinLength)
		if(nrow(CNVs) > 0)
		{
			ptm.tmp <- proc.time()
			CNVsRes <- FilterCNVs.V4(CNVs = CNVs, MinNumSNPs=MINNumSNPs, CNV, ID) # PathRawData = PathRawData,
			Res.tmp <- proc.time() - ptm.tmp
			#cat("Filter CNVs time: ", Res.tmp["elapsed"], "\n")
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
	df$CN <- as.numeric(df$CN)
		

	TimeRes <- proc.time() - ptm
	cat("Total time: ", TimeRes["elapsed"], "\n")
	return(df)
}
