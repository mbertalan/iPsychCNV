##' CheckAbnChr: Estimate abnormal chromosome.
##'
##' Specifically designed to reduce false positive CNVs and handle data from amplified DNA on dried blood spots.
##' @title CheckAbnChr
##' @param Path2RawFiles: Path for the Log R Ratio (LRR) and B Allele Frequency (BAF) files. Example: "/media/NeoScreen/NeSc_home/ILMN/iPSYCH/Version2".
##' @param Cores: Number of cores to run in parallel, default = 1. 
##' @param Pattern: Files pattern in the path. Example: "*.txt$".
##' @param Skip: Integer, the number of lines of the data file to be skipped before beginning to read data. Use if file has comments, default = 10.
##' @param NumFiles: Number of files to run. Example: numeric, 10 or character "All", default = All. 
##' @return Data frame with the estimate copy number for each chromosome.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' mockCNV <- MockData(N=5, Type="Blood", Cores=1)
##' cnvs <- CheckAbnChr(PathRawData=".", Cores=1, Pattern="^MockSample*", Skip=0)

CheckAbnChr <- function(Path2RawFiles="/media/NeoScreen/NeSc_home/ILMN/iPSYCH/Version2", Files=NA, Cores=1, Pattern="*.txt$", skip=10, NumFiles="All")
{
	library(iPsychCNV)
	library(parallel)
	
	suppressWarnings(if(is.na(Files))
	{
		Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=recursive)
	})
	
	if(NumFiles %in% "All"){ NumFiles <- length(Files) }
		
	Res2 <- mclapply(Files[1:NumFiles], mc.cores=Cores, mc.preschedule = FALSE, function(X)
	{
		File <- X
		SampleID <- tail(unlist(strsplit(X, "/")),n=1)
		cat(File, "\n")

		# Reading sample
		Sample <- ReadSample(File, skip=skip)
		LocalID <- tail(unlist(strsplit(File, "/")), n=1)
		#LocalID <- unique(Sample$Sample.ID)[1]

		Res <- sapply(unique(Sample$Chr), function(CHR)
		{
			Sample.Y <- subset(Sample, Chr %in% CHR)
			Sample.Y <- Sample.Y[order(Sample.Y$Position),]
		
			# Allele # deCODE has no Allele
			#Allele.Y <- apply(Sample.Y, 1, function(X){ X["Allele1"] %in% X["Allele2"] })
			#Allele.Y.Perc <- (sum(Allele.Y)/length(Allele.Y))*100

			# BAF
			Class.Y <- ClassNumbers(Sample.Y)
			LRRmean.Y <- mean(Sample.Y$Log.R.Ratio)

			chrY.Perc <- (sum(Sample.Y$B.Allele.Freq > 0.2 & Sample.Y$B.Allele.Freq < 0.8)/nrow(Sample.Y))*100

			# Classification
			res2 <- data.frame(CNVmean=LRRmean.Y)
			BAF.Y <- EvaluateMyBAF(Class.Y, res2)
		
			Class.Y[1:7] <- Class.Y[1:7]
	
			df <- data.frame(LocalID=LocalID, Chr=CHR, chr.Perc=chrY.Perc,  LRRmean=LRRmean.Y, SampleID=SampleID, NumChr=BAF.Y, stringsAsFactors=F)
			df <- cbind(df, Class.Y)
			return(df)
		})
		df <- MatrixOrList2df(Res)
		
		# EvalChr
		df$ChrEval <- NA
		if(df$NumChr[df$Chr %in% "X"] == 1 & df$NumChr[df$Chr %in% "Y"] == 1){ df$ChrEval <- "Male" }
		if(df$NumChr[df$Chr %in% "X"] == 2 & df$NumChr[df$Chr %in% "Y"] == 0){ df$ChrEval <- "Female" }
		if(df$NumChr[df$Chr %in% "X"] == 1 & df$NumChr[df$Chr %in% "Y"] == 0){ df$ChrEval <- "Turner syndrome" }
		if(df$NumChr[df$Chr %in% "X"] == 3 & df$NumChr[df$Chr %in% "Y"] == 0){ df$ChrEval <- "Triple-X syndrome" }		
		if(df$NumChr[df$Chr %in% "X"] == 2 & df$NumChr[df$Chr %in% "Y"] == 1){ df$ChrEval <- "Klinefelter syndrome" }
		if(df$NumChr[df$Chr %in% "X"] == 1 & df$NumChr[df$Chr %in% "Y"] == 2){ df$ChrEval <- "XYY syndrome" } 
		if(df$NumChr[df$Chr %in% "X"] == 2 & df$NumChr[df$Chr %in% "Y"] == 2){ df$ChrEval <- "XXYY syndrome" }
		
		if(df$NumChr[df$Chr %in% "12"] == 3){ df$ChrEval <- "Trisomia 12" }
		if(df$NumChr[df$Chr %in% "13"] == 3){ df$ChrEval <- "Trisomia 13" }
		if(df$NumChr[df$Chr %in% "16"] == 3){ df$ChrEval <- "Trisomia 16" }
		if(df$NumChr[df$Chr %in% "17"] == 3){ df$ChrEval <- "Trisomia 17" }
		if(df$NumChr[df$Chr %in% "18"] == 3){ df$ChrEval <- "Trisomia 18" }
		if(df$NumChr[df$Chr %in% "21"] == 3){ df$ChrEval <- "Trisomia 21" }
		if(df$NumChr[df$Chr %in% "22"] == 3){ df$ChrEval <- "Trisomia 22" }

		return(df)
	})

	df <- MatrixOrList2df(Res2)
	return(df)
}
