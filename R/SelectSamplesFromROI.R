##' SelectSamplesFromROI: Select region of interest based on the samples.
##'
##' Given hotspots regions and overlap cutoff it will return CNVs that belong to each hotspot.
##' @title SelectSamplesFromROI
##' @param DF: Data frame with predicted CNVs for each sample. 
##' @param Roi: Regions Of Interest or hotspots.
##' @param OverlapMin: Minimum overlap a CNV need to have with the hotspot to be selected, default = 0.8. 
##' @param OverlapMax: Maximum overlap accepted for a CNV and a hotspots, default = 1.2.  
##' @return Data frame with predicted CNVs that belong to each hotspot.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @export
##' @examples
##' cnvs <- SelectSamplesFromROI(DF=Hotspots, roi=roi, OverlapMin=0.8, OverlapMax=1.2)


SelectSamplesFromROI <- function(DF, roi, OverlapMin=0.8, OverlapMax=1.2, Cores=1)
{
	suppressPackageStartupMessages(library(parallel))

	tmp2 <- mclapply(1:nrow(roi), mc.cores=Cores, mc.preschedule = FALSE, function(i) 
	{
		X <- roi[i,]
		LengthRoi <- as.numeric(X["Length"])
		StartRoi <- as.numeric(X["Start"])
		StopRoi <- as.numeric(X["Stop"])
		ChrRoi <-gsub(" ", "", X["Chr"])
		#Locus <- as.character(X["locus"])
		ROI <- paste(ChrRoi, StartRoi, StopRoi, sep="_", collapse="")
		cat(paste(i,"/", nrow(roi), "Chr:", ChrRoi, "Start:", StartRoi, "Stop:", StopRoi, sep="\t", collapse="\t"), "\r") 

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
				
				# Calculate the overlap region, than the percentage with ref and with CNV.
				if(Start >= StartRoi & Stop >= StopRoi) # |-->--|---<
				{				     
					   OverlapRef <- ((StopRoi - Start)/LengthRoi)
					   OverlapCNV <- ((StopRoi - Start)/Length)
				}
				if(Start <= StartRoi & Stop <= StopRoi) # >---|----<----|
				{					
					OverlapRef <- ((Stop - StartRoi)/LengthRoi)
					OverlapCNV <- ((Stop - StartRoi)/Length)
				}
				if(Start >= StartRoi & Stop <= StopRoi) # |--->-----<---|
				{
					OverlapRef <- (Length/LengthRoi)
					OverlapCNV <- Length/Length
				}
				if(Start <= StartRoi & Stop >= StopRoi) # >---|-------|---<
				{
					OverlapRef <- (Length/LengthRoi)
					OverlapCNV <- Length/Length
				}
				df <- data.frame(OverlapRef=OverlapRef, OverlapCNV=OverlapCNV, ROI=ROI, stringsAsFactors=FALSE)
				return(df)
			})
			Overlap <- MatrixOrList2df(Overlap)
			tmp2 <- cbind(tmp, Overlap)
			#tmp2$locus <- rep(Locus, nrow(tmp2))
			tmp <- subset(tmp2, OverlapRef > OverlapMin & OverlapRef < OverlapMax & OverlapCNV > OverlapMin & OverlapCNV < OverlapMax)
			#cat(nrow(tmp), "\n")
			return(tmp)
		}
	})
	tmp3 <- MatrixOrList2df(tmp2)
	return(tmp3)
}
