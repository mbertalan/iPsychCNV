RunPennCNV <- function(PathRawData = "~/CNVs/MockData/PKU/Data", Pattern=".*Mock.*\\.tab$", Cores=20, Skip=0, Normalization=FALSE, PFB=NA, HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/")

{
	library(parallel)
	Files <- list.files(path=PathRawData, pattern=Pattern, full.names=TRUE, recursive=FALSE)
	# Re-writing file.
	mclapply(Files, mc.cores=Cores, mc.preschedule = FALSE, function(file) 
	{
		# for some odd reason penncnv needs a name before LRR and BAF.
		cat(file, "\n")
		tmp <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=F)
		colnames(tmp)[colnames(tmp) %in% "B.Allele.Freq"] <- "C B Allele Freq"
		colnames(tmp)[colnames(tmp) %in% "Log.R.Ratio"] <- "C Log R Ratio"
		colnames(tmp)[colnames(tmp) %in% "SNP.Name"] <- "Name"
		tmp <- tmp[, c("Name", "Position", "Chr", "C B Allele Freq", "C Log R Ratio")]
		Newfile <- paste(file, ".penncnv", sep="", collapse="")
		write.table(tmp, file=Newfile, quote=FALSE, row.names=FALSE, sep="\t")
	})
	Files <- list.files(path=PathRawData, pattern=".*\\.penncnv$", full.names=TRUE, recursive=FALSE)
	
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
	

	Res <- mclapply(Files, mc.cores=Cores, mc.preschedule = FALSE, function(X) 
	{
		cat(X, "\n")
		ID <- tail(unlist(strsplit(X, "/")),n=1)
		
		if(Normalization)
		{
			CNV <- NormalizeData(CNV, ExpectedMean=0, DF=NA, FALSE)
		}
		
		# Find CNVs
		Output <- paste(ID, ".penncnv.out", sep="", collapse="")
		Command <- paste(Path2PennCNV, "detect_cnv.pl -test -minsnp 28 --minlength 10 --confidence -hmm ", HMM, " -pfb ", PFB, " ", X, " -log logfile -out ", Output, sep="", collapse="")
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
		ID2 <- sapply(ID, function(X){ unlist(strsplit(X, ".penncnv."))[1]  })
		tmp$ID <- ID2 
		tmp$CNVID <- 1:nrow(tmp)
		df <- tmp
		
		# Getting CNV mean from sample.
		Sample <- read.table(X, sep="\t", header=TRUE, stringsAsFactors=F)
		CNVmean <- apply(df[,c("Start", "Stop", "Chr")], 1, function(X)
		{
		 	Start <- X["Start"]
		 	Stop <- X["Stop"]
		 	CHR <- X["Chr"]
		 	subSample <- subset(Sample, Chr %in% CHR) 
		 	LRR <- subSample[,grep("Log.R.Ratio", colnames(subSample))]
		 	CNVmean <- mean(LRR[Start:Stop])
		 	return(CNVmean)
		})
		df$CNVmean <- CNVmean
		
		return(df)
	})
	df <- MatrixOrList2df(Res)
	return(df)
}
