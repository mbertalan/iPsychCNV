FilterFromCNVs <- function(CNVs, PathRawData, MinNumSNPs=10, Source="iPsychCNV", Skip=10, Cores)
{
	suppressPackageStartupMessages(library(parallel))
	
	tmp <- mclapply(unique(CNVs$ID), mc.cores=Cores, mc.preschedule = FALSE, function(IDs) 
	{
		subCNVs <- subset(CNVs, ID %in% IDs)
		#RawFile <- paste(PathRawData, "/", IDs, sep="", collapse="")
		RawFile <- IDs
		cat(IDs, "\n")
		Sample <- ReadSample(RawFile, skip=Skip, LCR=FALSE)

		Res <- apply(subCNVs, 1, function(X)
		{
			#cat(X, "\n")
			StartM <- as.numeric(X["Start"])
			StopM <- as.numeric(X["Stop"])
			ChrM <- X["Chr"]
			ChrM <- gsub(" ","", ChrM)

			subSample <- subset(Sample, Chr %in% ChrM)
			#cat(nrow(subSample), "\n")
			subSample <- subSample[with(subSample, order(Position)),]
			#NewName <- paste("Test", ChrM, StartM, StopM, ".RData", sep="", collapse="")
			#save(subSample, Sample, file=NewName)
			StartIndx <- which(subSample$Position == StartM)[1] 
			StopIndx <- which(subSample$Position == StopM)[1]
			Vector <- c(StartIndx, StopIndx, X)
			cat(ChrM, StartM, StopM, StartIndx, StopIndx, RawFile, "\n")
			names(Vector)[1:2] <- c("StartIndx", "StopIndx")
			return(Vector)
		})
		tmp <- as.data.frame(t(Res), stringsAsFactors=F)

		tmp2 <- FilterCNVs.V4(tmp, MinNumSNPs=MinNumSNPs, CNV=Sample, ID=IDs)
		return(tmp2)
	})
	tmp2 <- MatrixOrList2df(tmp)
	df <- tmp2[, !colnames(tmp2) %in% ".id"]
	df$Source <- rep(Source, nrow(df))
	df$CN <- df$Class
	df$CN[df$CN %in% "Del"] <- "1"
	df$CN[df$CN %in% "Normal"] <- "2"
	df$CN[df$CN %in% "Dup"] <- "3"
	df$CN[df$CN %in% "DoubleDel"] <- "0"
	df$CN[df$CN %in% "DoubleDup"] <- "4"
	df$CN <- as.numeric(df$CN)
	return(df)
}
