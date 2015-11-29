##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @return Data frame with CNVs predicted.
##' @author Marcelo Bertalan
##' @export
RunPennCNV <- function(PathRawData = "~/CNVs/MockData/PKU/Data", MINNumSNPs=20, Pattern=".*Mock.*\\.tab$", Cores=20, Skip=0, Normalization=FALSE, PFB="NO", HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/",  penalty=60, Quantile=TRUE, QSpline=TRUE, sd=0.15, PennCNVFormat=FALSE, RemoveTmpfiles=FALSE)
{
	library(parallel)
	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=FALSE)

	# Re-writing file.
	if(!PennCNVFormat) # If not in PennCNV format, than make files to match PennCNV format.
	{
		mclapply(Files, mc.cores=Cores, mc.preschedule = FALSE, function(file) 
		{
			# penncnv needs a name before LRR and BAF.
			cat(file, "\n")
			tmp <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=F, skip=Skip)
			tmp <- subset(tmp, !is.na(tmp$B.Allele.Freq))
			tmp <- subset(tmp, !is.na(tmp$Log.R.Ratio))
		
			if(Normalization)
			{
				tmp <- NormalizeData(tmp,ExpectedMean=0, penalty=penalty, Quantile=Quantile, QSpline=QSpline, sd=sd)
			}
		
			colnames(tmp)[colnames(tmp) %in% "B.Allele.Freq"] <- "C B Allele Freq"
			colnames(tmp)[colnames(tmp) %in% "Log.R.Ratio"] <- "C Log R Ratio"
			colnames(tmp)[colnames(tmp) %in% "SNP.Name"] <- "Name"
			tmp <- tmp[, c("Name", "Position", "Chr", "C B Allele Freq", "C Log R Ratio")]
			Newfile <- paste(file, ".penncnv", sep="", collapse="")
			write.table(tmp, file=Newfile, quote=FALSE, row.names=FALSE, sep="\t")
		})
		Files <- list.files(path=PathRawData, pattern=".*\\.penncnv$", full.names=TRUE, recursive=FALSE)
	}
	
	# Creating Chip info file for PennCNV
	cat(Files[1], "\n")
	ChipInfo <- read.table(Files[1], sep="\t", header=TRUE, stringsAsFactors=F)
	ChipInfo <- ChipInfo[, c("Name", "Chr", "Position")]
	write.table(ChipInfo, file="SNP.Position.tab", quote=FALSE, row.names=FALSE, sep="\t")

	write.table(Files, "Mock.List.Files.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
		
	if(is.na(PFB))
	{
		Command <- paste(Path2PennCNV, "compile_pfb.pl -listfile ./Mock.List.Files.txt -output Mock.pfb", sep="", collapse="")
		cat(Command, "\n")
		system(Command)
		PFB <- "Mock.pfb"
	}
	else if(PFB %in% "NO")
	{
		ChipInfo$PFB <- 0.5
		write.table(ChipInfo, file="Mock.pfb", quote=FALSE, row.names=FALSE, sep="\t")
		PFB <- "Mock.pfb"
	}
		
		
	

	Res <- mclapply(Files, mc.cores=Cores, mc.preschedule = FALSE, function(X) 
	{
		cat(X, "\n")
		ID <- tail(unlist(strsplit(X, "/")),n=1)
		
		# Find CNVs
		Output <- paste(ID, ".penncnv.out", sep="", collapse="")
		Command <- paste(Path2PennCNV, "detect_cnv.pl -test -minsnp ", MINNumSNPs, " --minlength 10 --confidence -hmm ", HMM, " -pfb ", PFB, " ", X, " -log logfile -out ", Output, sep="", collapse="")
		cat(Command, "\n")
		system(Command)
		
		# Merge CNVs
		OutputMerged <- paste(ID, ".penncnv.out.merged", sep="", collapse="")
		Command <- paste(Path2PennCNV, "clean_cnv.pl combineseg ", Output, " --signalfile SNP.Position.tab --fraction 0.2 --bp --output ", OutputMerged, sep="", collapse="")
		cat(Command, "\n")
		system(Command)
		Penn2Tab <- system.file("exec/Penn2Tab.pl",package="iPsychCNV")
		
		# Change format
		OutputMergedTab <- paste(ID, ".penncnv.out.merged.tab", sep="", collapse="")
		Command <- paste(Penn2Tab, " < ", OutputMerged, " > ", OutputMergedTab, sep="", collapse="")
		cat(Command, "\n")
		system(Command)
		
		# Reading the result and adding ID.
		tmp <- read.table(OutputMergedTab, sep="\t", header=TRUE, stringsAsFactors=FALSE)
		ID <- sapply(tmp$File, function(X){ ID <- tail(unlist(strsplit(X, "/")),n=1)  })
		ID2 <- sapply(ID, function(X){ unlist(strsplit(X, ".penncnv"))[1]  })
		tmp$ID <- ID2 
		tmp$CNVID <- 1:nrow(tmp)
		df <- tmp
		
		# Getting CNV mean from sample.
		Sample <- read.table(X, sep="\t", header=TRUE, stringsAsFactors=F)
		CNVmean <- apply(df[,c("Start", "Stop", "Chr")], 1, function(X)
		{
		 	Start <- as.numeric(X["Start"])
		 	Stop <- as.numeric(X["Stop"])
		 	CHR <- X["Chr"]
		 	subSample <- subset(Sample, Chr %in% CHR & Position >= Start & Position <= Stop) 
		 	LRR <- subSample[,grep("Log.R.Ratio", colnames(subSample))]
		 	CNVmean <- mean(LRR)
		 	return(CNVmean)
		})
		df$CNVmean <- CNVmean
		
		return(df)
	})
	
	df <- MatrixOrList2df(Res)
	df$Source <- "PennCNV"
	
	# Removing tmp files
	if(RemoveTmpfiles)
	{
		RemoveTmpFiles <- paste(PathRawData, "/*penncnv*", sep="", collapse="")
		unlink(x=RemoveTmpFiles)
	}
	
	return(df)
}
