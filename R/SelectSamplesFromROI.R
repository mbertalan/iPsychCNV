##' SelectSamplesFromROI
##'
##' Given hotspots regions and overlap cutoff it will return CNVs that belong to each hotspot.
##' @title SelectSamplesFromTOI
##' @return Data frame with CNVs predicted that belong to each hotspot.
##' @param DF: Data frame with CNVs predicted for each sample. 
##' @param roi: Regions of Interest or hotspots.
##' @param OverlapMin: Minimum overlap a CNV need to have with the hotspot to be selected. 
##' @param OverlapMax: Maximum overlap accepted for a CNV and a hotspots.  
##' @author Marcelo Bertalan
##' @export
##' @examples
##' cnvs <- SelectSamplesFromROI(DF=Hotspots, roi=roi, OverlapMin=0.8, OverlapMax=1.2)


SelectSamplesFromROI <- function(DF, roi, OverlapMin=0.8, OverlapMax=1.2, Cores=1)
{
	suppressPackageStartupMessages(library(parallel))
	
	tmp2 <- mclapply(1:nrow(roi), mc.cores=Cores, mc.preschedule = FALSE, function(i) 
	{
		X <- roi[i,]
		print(X, sep="\r")
		LengthRoi <- as.numeric(X["Length"])
		StartRoi <- as.numeric(X["Start"])
		StopRoi <- as.numeric(X["Stop"])
		ChrRoi <-gsub(" ", "", X["Chr"])
		Locus <- as.character(X["locus"])
		tmp <- subset(DF, Chr %in% ChrRoi & Start < StopRoi & Stop > StartRoi)  # & ((Length/LengthRoi) >  Overlap)
		if(nrow(tmp) > 0)
		{
			#cat(nrow(tmp), "\n")
			#Overlap
			Overlap <- apply(tmp, 1, function(Y)
			{
				Length <- as.numeric(Y["Length"])
				Start <- as.numeric(Y["Start"])
				Stop <- as.numeric(Y["Stop"])
				#cat(Start, Stop, StartRoi, StopRoi, "\n")
				if(Start >= StartRoi & Stop >= StopRoi) # |-->--|---<
				{				     
					   Overlaplength <- ((StopRoi - Start)/LengthRoi)
					   OverlapAll <- Length/LengthRoi
				}
				if(Start <= StartRoi & Stop <= StopRoi) # >---|----<----|
				{					
					Overlaplength <- ((Stop - StartRoi)/LengthRoi)
					OverlapAll <- Length/LengthRoi
				}
				if(Start >= StartRoi & Stop <= StopRoi) # |--->-----<---|
				{
					Overlaplength <- (Length/LengthRoi)
					OverlapAll <- Length/LengthRoi
				}
				if(Start <= StartRoi & Stop >= StopRoi) # >---|-------|---<
				{
					Overlaplength <- (Length/LengthRoi)
					OverlapAll <- Length/LengthRoi
				}
				df <- data.frame(Overlaplength=Overlaplength, OverlapAll=OverlapAll)
				return(df)
			})
			Overlap <- MatrixOrList2df(Overlap)
			tmp2 <- cbind(tmp, Overlap)
			tmp2$Locus <- rep(Locus, nrow(tmp2))
			tmp <- subset(tmp2, OverlapAll > OverlapMin & OverlapAll < OverlapMax & Overlaplength > OverlapMin)
			#cat(nrow(tmp), "\n")
			return(tmp)
		}
	})
	tmp3 <- MatrixOrList2df(tmp2)
	return(tmp3)
}
