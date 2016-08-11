##' GetMockValues: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iGetMockValues
##' @param Type: Unkown, default = Blood. 
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

GetMockValues <- function(Type="Blood")
{
	# Telomere noise
	TelomereNoiseSize <- seq(from=100, to=1300, by=100)
	TelomereNoiseEffect <- seq(from=0.5, to=-0.7, by=-0.1)

	# BadSNPs      1     2    3    4    5    6   7     8    9   10   11   12   13   14   15   16   17   18   19   20   21   22
	BadSNPs <- c(0.02, 0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.001,0.02,0.1,0.1,0.15,0.001,0.25,0.05,0.05,0.15)
	names(BadSNPs) <- 1:22
	BadSNPIntensity <- seq(from=0.2, to=-4, by=-0.1)
	BadSNPIntensityProb <- seq(from=0.53, to=0.11, by=-0.01)

	# BAF INFO
	BAFs <- seq(from=0, to=1, by=0.01) # 101
	BAF_Basic <- rep(0.00001, 101)
	names(BAF_Basic) <- BAFs

	BAF_Normal <- BAF_Basic
	BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + 1
	BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + 0.05
	
	BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + 0.05
	BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + 1
	
	BAF_Normal[c(47:53)] <- BAF_Normal[c(47:53)] + 0.05
	BAF_Normal[c(50:51)] <- BAF_Normal[c(50:51)] + 0.1

	# BAF Del prob
	BAF_Del <- BAF_Basic
	BAF_Del[1:2] <- BAF_Del[1:2] + 1
	BAF_Del[100:101] <- BAF_Del[100:101] + 1
	
	# BAF Dup prob
	BAF_Dup <-  BAF_Basic
	BAF_Dup[1:2] <- BAF_Dup[1:2] + 1
	BAF_Dup[100:101] <- BAF_Dup[100:101] + 1
	BAF_Dup[30:35] <- BAF_Dup[30:35] + 0.05 
	BAF_Dup[32:33] <- BAF_Dup[32:33] + 0.1 
	BAF_Dup[65:70] <- BAF_Dup[65:70] + 0.05
	BAF_Dup[67:68] <- BAF_Dup[67:68] + 0.1 	

	# BAF CN=0
	BAF_CN0 <- BAF_Basic
	
	# BAF CN=4
	BAF_CN4 <- BAF_Basic
	
	BAF_CN4[1:2] <- BAF_CN4[1:2] + 1
	BAF_CN4[100:101] <- BAF_CN4[100:101] + 1
	
	BAF_CN4[21:29] <- BAF_CN4[21:29] + 0.05 
	BAF_CN4[23:27] <- BAF_CN4[23:27] + 0.1
	
	BAF_Dup[71:79] <- BAF_Dup[71:79] + 0.05
	BAF_Dup[73:77] <- BAF_Dup[73:77] + 0.1
	
	BAF_CN4[c(48:53)] <- BAF_CN4[c(48:53)] + 0.05
	BAF_CN4[c(50:51)] <- BAF_CN4[c(50:51)] + 0.1
	

	if(Type %in% "Blood")
	{
		# CNV Info
		Del <- seq(from=-0.3, to=-0.5, by=-0.05)
		Dup <- seq(from=0.3, to=0.5, by=0.05)
		ChrSD <- seq(from=0.1, to=0.5, by=0.025)
		ChrSDProb <- c(0.4, 0.9, 0.8, 0.4, 0.4, 0.2, 0.1, 0.01,0.01, 0.01,0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
		ChrMean <- ChrMean_Blood_800_Density
		ChrMeanProb <- ChrMeanProb_Blood_800_Density
	}
	else # PKU
	{
		# CNV Info
		Del <- seq(from=-0.15, to=-0.5, by=-0.05)
		Dup <- seq(from=0.15, to=0.5, by=0.05)
		ChrSD <- seq(from=0.1, to=0.5, by=0.025)
		ChrSDProb <- c(0.2, 0.6, 0.9, 0.6, 0.4,  0.2, 0.25, 0.2, 0.15, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01)
		ChrMean <- ChrMean_PKU_Wave10_Density
		ChrMeanProb <- ChrMeanProb_PKU_Wave10_Density
	}
	tmp <- list(Del=Del, Dup=Dup, BAF_CN0=BAF_CN0, BAF_CN4=BAF_CN4, ChrMean=ChrMean, ChrSD=ChrSD, ChrSDProb=ChrSDProb, TelomereNoiseSize=TelomereNoiseSize, TelomereNoiseEffect=TelomereNoiseEffect, BAF_Normal=BAF_Normal, BAF_Del=BAF_Del, BAF_Dup=BAF_Dup, BadSNPs=BadSNPs, BadSNPIntensity=BadSNPIntensity, BadSNPIntensityProb=BadSNPIntensityProb, BAFs=BAFs, ChrMeanProb=ChrMeanProb)
	return(tmp)
}
