##' DefineStartAndStop: Define Start and Stop for CNV. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title DefineStartAndStop
##' @return Data frame with Start and Stop.
##' @author Marcelo Bertalan
##' @export

DefineStartAndStop <- function(indx, subCNV, MinNumSNPs, CHR, ID, CPT.Res)
{
	library(pastecs)
	
	Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
	Info <- t(Info)
	Df <- as.data.frame(Info)
	CNVMean <- apply(Df, 1, function(X){ mean(CPT.Res@data.set[X[1]:X[2]]) } )
	Df$CNVMean <- CNVMean
	names(Df) <- c("Start", "Stop", "CNVMean")


	Df.fixpos <- apply(Df, 1, function(X){ range(which(abs(CPT.Res@data.set[X[1]:X[2]]) > 0.1))+(X[1]-1) }) # 0.05 or 0.1 ? Ok if edges are low values
	DF <- as.data.frame(t(Df.fixpos))
	names(DF) <- c("Start", "Stop") # Df with indx
	DF <- DF[is.finite(DF[,1]),]
	if(nrow(DF) > 0)
	{
		DF$NumSNPs <- as.numeric(DF$Stop) - as.numeric(DF$Start)
		DF <- subset(DF, NumSNPs > MinNumSNPs)
	}

	# Check for pits
	if(nrow(DF) > 0)
	{
		NewDF <- apply(DF, 1, function(X)
		{
			tp <- turnpoints(abs(CPT.Res@data.set[X[1]:X[2]]))
			Pits <- which(tp$pits)
			Peaks <- which(tp$peaks)
			if(length(Pits) > 0 && length(Peaks) > 0) # Only change position if pits is before peak. Very common have a pit in the middle of the CNV.
			{ 
				if(Pits[1] < Peaks[1]){ NewStart <- X[1] + Pits[1] }else{ NewStart <- X[1] }
				if(tail(Pits, n=1) > tail(Peaks, n=1)){ NewStop <- (X[1]-1) + tail(Pits, n=1) }else{ NewStop <- X[2] }
			}
			else # No need for cut
			{ 
				NewStart <- X[1] 
				NewStop <- X[2] 
			}
			NewDf <- data.frame(Start=NewStart, Stop=NewStop)
			return(NewDf)
		})
	
		if(length(NewDF) > 0)
		{
			NewDF <- MatrixOrList2df(NewDF)
		}

		if(nrow(NewDF) > 0)
		{
			Df.fixpos <- apply(NewDF, 1, function(X){ range(which(abs(CPT.Res@data.set[X[1]:X[2]]) > 0.1))+X[1] }) 
			DF2 <- as.data.frame(t(Df.fixpos))
			DF2 <- DF2[,1:2]
			names(DF2) <- c("Start", "Stop") # Df with indx
	
			DF3 <- Plot.CNV.Info(MinNumSNPs, DF2, subCNV, ID) # Min number of SNPs for CNV
			
			return(DF3)
		}
	}	
}
