SelectSamplesFromROI <- function(DF, roi, OverlapMin, OverlapMax)
{
	tmp2 <- apply(roi, 1, function(X)
	{
		cat(X, "\n")
		LengthRoi <- as.numeric(X["Length"])
		StartRoi <- as.numeric(X["Start"])
		StopRoi <- as.numeric(X["Stop"])
		ChrRoi <- as.character(X["Chr"])
		tmp <- subset(DF, Chr %in% ChrRoi & Start < StopRoi & Stop > StartRoi)  # & ((Length/LengthRoi) >  Overlap)
		if(nrow(tmp) > 0)
		{
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
				}
				if(Start <= StartRoi & Stop <= StopRoi) # >---|----<----|
				{					
					Overlaplength <- ((Stop - StartRoi)/LengthRoi)
				}
				if(Start >= StartRoi & Stop <= StopRoi) # |--->-----<---|
				{
					Overlaplength <- (Length/LengthRoi)
				}
				if(Start <= StartRoi & Stop >= StopRoi) # >---|-------|---<
				{
					Overlaplength <- (Length/LengthRoi)
				}
				return(Overlaplength)
			})
			tmp$Overlap <- Overlap
			tmp <- subset(tmp, Overlap > OverlapMin & Overlap < OverlapMax)
			return(tmp)
		}
	})
	tmp3 <- MatrixOrList2df(tmp2)
	return(tmp3)
}
