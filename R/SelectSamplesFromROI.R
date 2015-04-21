SelectSamplesFromROI <- function(DF, roi, OverlapMin, OverlapMax)
{
	tmp2 <- apply(roi, 1, function(X)
	{
		cat(X, "\n")
		LengthRoi <- as.numeric(X["Length"])
		StartRoi <- as.numeric(X["Start"])
		StopRoi <- as.numeric(X["Stop"])
		ChrRoi <- as.character(X["Chr"])
		Locus <- as.character(X["locus"])
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
			tmp <- subset(tmp2, OverlapAll > OverlapMin & OverlapAll < OverlapMax)
			return(tmp)
		}
	})
	tmp3 <- MatrixOrList2df(tmp2)
	return(tmp3)
}
