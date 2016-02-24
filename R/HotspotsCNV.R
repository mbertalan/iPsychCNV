##' HotspotsCNV: Identify Copy Number Variation (CNV) hotspots. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title HotspotsCNV
##' @param Df: Dataframe with CNV predction with chromosome (Chr.), start, and stop position.
##' @param Freq: Minimum number of CNVs to be considered a hotspot, default = 1.
##' @param OverlapCutoff: Minimum overlap among CNVs to be considered the same CNV region, default = 0.7.
##' @param Cores: Numeric, Number of cores used, default = 1. 
##' @param OverlapMin: Minimum overlap with hotspot to be selected for counting, default = 0.9.
##' @param OverlapMax: Maximum overlap with hotspot to be selected for counting, default = 1.1.
##' @return CNV hotspots.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' MockDataCNVs <- MockData(N=100, Type="PKU", Cores=20)
##' iPsych.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=20, Cores=20, Pattern="^MockSample", MinLength=10, Skip=0)
##' iPsych.Pred.hotspots <- HotspotsCNV(iPsych.Pred, Freq=2, OverlapCutoff=0.9, Cores=1)

HotspotsCNV <- function(df, Freq=1, OverlapCutoff=0.7, Cores=1, OverlapMin=0.9,  OverlapMax=1.1)
{
	library(plyr)
	library(parallel)

	if(length(df$Source[1]) > 0){ Source <- df$Source[1] }else{ stop("Data frame does not have Source\n") }
	OriginalDF <- df

	tmp <- apply(df, 1, function(X)
	{
		Start1 <- as.numeric(X["Start"])
		Stop1 <- as.numeric(X["Stop"])
		Chr1 <- gsub(" ", "", X["Chr"])
		Length1 <- as.numeric(X["Length"])
		CNV_Pos <- paste(Chr1, Start1, Stop1, Length1, sep="_", collapse="") 
		return(CNV_Pos)
	})

	tmp2 <- table(tmp)

	# selecting by freq
	tmp3 <- tmp2[tmp2 > Freq]
	tmp4 <- as.vector(tmp3)
	names(tmp4) <- names(tmp3)

	# selecting by size
	tmp <- sapply(1:length(tmp4), function(X)
	{
		cat(X, "\r")
		df1 <- unlist(strsplit(names(tmp4)[X], "_"))
		df2 <- c(df1, tmp4[X], names(tmp4)[X])
		names(df2) <- c("Chr", "Start", "Stop", "Length", "Count", "NewName")
		return(df2)
	})
	tmp4 <- as.data.frame(t(tmp), stringsAsFactors=F)
	tmp4$Source <- rep(Source, nrow(tmp4))
	
	cat("Running CompressCNVs", nrow(tmp4), "\n")
	tmp5 <- CompressCNVs(tmp4, OverlapCutoff, Cores) # Problem
	cat("Re-running CompressCNVs", nrow(tmp5), "\n")
	NumOfCNVs <- 0
	while(nrow(tmp5) != NumOfCNVs)
	{
		NumOfCNVs <- nrow(tmp5)
		tmp5 <- CompressCNVs(tmp5, OverlapCutoff, Cores)
		if(nrow(tmp5) == NumOfCNVs)
		{
			cat("CompressCNVs done", nrow(tmp5), NumOfCNVs,  "\n")
		}
		else
		{
			cat("Re-running CompressCNVs", nrow(tmp5), NumOfCNVs,  "\n")
		}
	}
	#save(tmp5, file="tmp5.RData")
	# Count CNVs in compressed CNV regions	
	cat("Counting CNVs\n")

	#save(tmp5, file="tmp5.RData")
	CNV_Count <- sapply(1:nrow(tmp5), function(X)
	{
		tmp <- SelectSamplesFromROI(DF=OriginalDF, roi=tmp5[X,], OverlapMin=OverlapMin,  OverlapMax=OverlapMax)
		Counts <- data.frame("CN0"=0, "CN1"=0, "CN2"=0, "CN3"=0, "CN4"=0, "CN5"=0)
		tmp2 <- as.data.frame(table(tmp$CN), stringsAsFactors=F)
		tmp3 <- tmp2$Freq
		names(tmp3) <- tmp2$Var1
		names(tmp3) <- sapply(names(tmp3), function(X){ paste("CN", X, sep="", collapse="") })
		Counts[names(tmp3)] <- tmp3
		return(Counts)
	})

	CNV_Count2 <- MatrixOrList2df(CNV_Count)
	#save(CNV_Count2, file="CNV_Count.RData")
	#save(tmp5, file="tmp5.RData")
	tmp5 <- GetLocus(tmp5)	
	tmp5$ID <- rep("ROI", nrow(tmp5))	
	tmp5$Class <- rep("ROI", nrow(tmp5))		
	tmp6 <- cbind(tmp5, CNV_Count2)
	tmp6$CNTotal <- tmp6$CN0 + tmp6$CN1 + tmp6$CN3 + tmp6$CN4
	return(tmp6)
}
