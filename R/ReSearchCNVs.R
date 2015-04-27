##' ReSearching: re-search all CNVs found in all samples.
##'
##' In noisy data the location of CNVs by Log.R.Ratio may fail. Using re-searching is possible to check if there are good CNVs in other samples. This will improve th eestimate frequency of a CNV in a population.
##' @title ReSearching.
##' @return Returns a data frame with CNVs.
##' @author Marcelo Bertalan
##' @export

ReSearchCNVs <- function(Files, Cores=1, Hotspots=roi, MINNumSNPs=20, chr=NA, Skip=10, LCR=TRUE, PFB=NULL, NormQspline=FALSE, Quantile=TRUE, DFspline=NA)
{
	tmp <- mclapply(Files, mc.cores=Cores, function(RawFile) 
	{
		cat("Re-Running:\t", RawFile, "\n")
		ID <- tail(unlist(strsplit(RawFile, "/")),n=1)
		Sample <- ReadSample(RawFile, skip=Skip, LCR=LCR, PFB=PFB, chr=chr)
		Sample <- NormalizeData(Sample, ExpectedMean=0, DF=DFspline, NormQspline=NormQspline, Quantile=Quantile)
		Hotspots <- Hotspots[,!colnames(Hotspots) %in% "ID"]

		Res <- apply(Hotspots, 1, function(X)
		{
			StartM <- as.numeric(X["Start"])
			StopM <- as.numeric(X["Stop"])
			ChrM <- X["Chr"]
			ChrM <- gsub(" ","", ChrM)

			subSample <- subset(Sample, Chr %in% ChrM)
			subSample <- subSample[with(subSample, order(Position)),]
			StartIndx <- which(subSample$Position == StartM)[1] 
			StopIndx <- which(subSample$Position == StopM)[1]
			NumSNPs <- StopIndx - StartIndx
			Vector <- c(StartIndx, StopIndx, NumSNPs, ID, X)
			names(Vector)[1:4] <- c("StartIndx", "StopIndx", "NumSNPs", "ID")
			return(Vector)
		})
		tmp <- as.data.frame(t(Res), stringsAsFactors=F)
		tmp2 <- FilterCNVs.V4(tmp, MinNumSNPs=MinNumSNPs, CNV=Sample, ID=ID)
	})
	df <- MatrixOrList2df(tmp)
	df$CN <- df$Class
	df$CN[df$CN %in% "Del"] <- "1"
	df$CN[df$CN %in% "Normal"] <- "2"
	df$CN[df$CN %in% "Dup"] <- "3"
	df$CN[df$CN %in% "DoubleDel"] <- "0"
	df$CN[df$CN %in% "DoubleDup"] <- "4"
	df$CN <- as.numeric(df$CN)
	return(df)
}
