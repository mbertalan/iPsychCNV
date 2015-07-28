##' CheckAbnChr
##'
##' @title CheckAbnChr
##' @return Abnormal chromosomes
##' @author Marcelo Bertalan
##' @export

CheckAbnChr <- function(Path2RawFiles="/media/NeoScreen/NeSc_home/ILMN/iPSYCH/Version2", Cores=40, Pattern="*.txt$", skip=10, NumFiles="All")
{
	library(iPsychCNV)
	library(parallel)
	Files <- list.files(Path2RawFiles, pattern=Pattern, recursive=TRUE) 
	Files <- Files[!is.na(Files)]

	if(NumFiles %in% "All"){ NumFiles <- length(Files) }
		
	Res2 <- mclapply(Files[1:NumFiles], mc.cores=Cores, mc.preschedule = FALSE, function(X)
	{
		SampleID <- X
		File <- paste(Path2RawFiles, "/", X, sep="", collapse="")
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
			LRRmean.Y <- median(Sample.Y$Log.R.Ratio)

			chrY.Perc <- (sum(Sample.Y$B.Allele.Freq > 0.1 & Sample.Y$B.Allele.Freq < 0.9)/nrow(Sample.Y))*100

			# Classification
			BAF.Y <- MyBAF(Class.Y, LRRmean.Y, chrY.Perc)
		
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
