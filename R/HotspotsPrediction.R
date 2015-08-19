##' HotspotsPrediction: Count number of true CNVs in hotspots. 
##'
##' @title HotspotsPrediction
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

HotspotsPrediction <- function(hotspots, df)
{ 	
	library(iPsychCNV)	
	Res <- apply(hotspots[,c("Start", "Stop", "Chr", "NumSNPs", "CN", "CNVmean", "sd")], 1, function(X)
	{
		StartH <- as.numeric(X[1]);
		StopH <- as.numeric(X[2]);
		ChrH <- X[3];
		NumSNPsH <- as.numeric(X[4]);
		CNH <- as.numeric(X[5]);
		CNVmeanH <- as.numeric(X[6]);
		sdH <- as.numeric(X[7]);
		ID <- paste(ChrH, StartH, StopH, sep=".", collapse="")

		tmp <- subset(df, Start < StopH & Stop > StartH & Chr %in% ChrH ) # & CN == CNH
		tmp$overlap <- tmp$NumSNPs/NumSNPsH
		tmp2 <- subset(tmp, overlap > 0.7)

		df <- data.frame(CNVs=nrow(tmp2), CN=CNH, NumSNPs=NumSNPsH, ID=ID, CNVmean=CNVmeanH, sd=sdH, stringsAsFactors=F)
		return(df)
	})
	tmp2 <- MatrixOrList2df(Res)
	return(tmp2)
}
