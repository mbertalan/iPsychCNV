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
			StartM <- as.numeric(X["Start"])
			StopM <- as.numeric(X["Stop"])
			ChrM <- X["Chr"]
			ChrM <- gsub(" ","", ChrM)

			subSample <- subset(Sample, Chr %in% ChrM)
			subSample <- subSample[with(subSample, order(Position)),]
			
			StartIndx <- which(subSample$Position == StartM)[1] 
			StopIndx <- which(subSample$Position == StopM)[1]
			df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx)
			cat(ChrM, StartM, StopM, StartIndx, StopIndx, RawFile, "\n")
			return(df)
		})
		df <- MatrixOrList2df(tmp=Res)
		tmp <- cbind(subCNVs, df)

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
