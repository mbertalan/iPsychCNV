GetMockValues <- function(Type="Blood", BAF_LOH=FALSE)
{
	# BAF Basic
	BAFs <- seq(from=0, to=1, by=0.01) # 101
	BAF_Basic <- rep(0.01, 101)
	names(BAF_Basic) <- BAFs
	
	# Telomere noise
	TelomereNoiseSize <- seq(from=100, to=1300, by=100)
	TelomereNoiseEffect <- seq(from=0.5, to=-0.7, by=-0.1)

	# BadSNPs      1     2    3    4    5    6   7     8    9   10   11   12   13   14   15   16   17   18   19   20   21   22
	BadSNPs <- c(0.02, 0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.001,0.02,0.1,0.1,0.15,0.001,0.25,0.05,0.05,0.15)
	names(BadSNPs) <- 1:22
	BadSNPIntensity <- seq(from=0.2, to=-4, by=-0.1)
	BadSNPIntensityProb <- seq(from=0.53, to=0.11, by=-0.01)

	# From BAF basic to Normal.
	if(BAF_LOH)
	{
		BAF_Basic[10:90] <- seq(from=0.18, to=0.22, by=0.0005)
	
		# BAF normal prob	
		BAF_Normal <- BAF_Basic
		BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + 0.2
		BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + 0.1
		BAF_Normal[c(8:9)] <- BAF_Normal[c(8:9)] - 0.02

		BAF_Normal[c(80:85)] <- BAF_Normal[c(80:85)] + 0.01
		BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + 0.20
		BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + 0.5	
	}
	else
	{		
		if(Type %in% "Blood")
		{
			BAF_Normal <- BAF_Basic
			BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + 0.52
			BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + 0.32
			BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + 0.8
			BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + 0.35
			BAF_Normal[c(45:55)] <- BAF_Normal[c(45:55)] + 0.05
			BAF_Normal[c(47:53)] <- BAF_Normal[c(47:53)] + 0.05
			BAF_Normal[c(50:51)] <- BAF_Normal[c(50:51)] + 0.1
		}
		else # PKU
		{
			BAF_Normal <- BAF_Basic + 0.02
			BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + 0.6
			BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + 0.3
			BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + 0.9
			BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + 0.55
			BAF_Normal[c(45:55)] <- BAF_Normal[c(45:55)] + 0.03
			BAF_Normal[c(47:53)] <- BAF_Normal[c(47:53)] + 0.05
			BAF_Normal[c(50:51)] <- BAF_Normal[c(50:51)] + 0.05
		}
	}

	if(Type %in% "Blood")
	{
		# CNV Info
		Del <- seq(from=-0.3, to=-0.5, by=-0.05)
		Dup <- seq(from=0.3, to=0.5, by=0.05)
		ChrMean <- seq(from=0.1, to=-0.1, by=-0.01)
		ChrMeanProb <- c(0.01, 0.01,0.01,0.01,0.02,0.02,0.1,0.3,0.3,0.15,0.1,0.02,0.2,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
		ChrSD <- seq(from=0.1, to=0.5, by=0.025)
		ChrSDProb <- c(0.2, 0.6, 0.7, 0.9, 0.7, 0.2, 0.1, 0.01,0.01, 0.01,0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
	
		# BAF Del prob
		BAF_Del <- BAF_Basic
		BAF_Del[1:2] <- BAF_Del[1:2] + 1
		BAF_Del[100:101] <- BAF_Del[100:101] + 1
		
		# BAF Dup prob
		BAF_Dup <-  BAF_Basic
		BAF_Dup[1:2] <- BAF_Dup[1:2] + 0.35
		BAF_Dup[100:101] <- BAF_Dup[100:101] + 0.35
		BAF_Dup[25:35] <- BAF_Dup[25:35] + 0.2 # 0.25 0.30 0.35, 
		BAF_Dup[70:80] <- BAF_Dup[70:80] + 0.2 # 0.7 0.75 0.8	
	}
	else
	{
		# CNV Info
		Del <- seq(from=-0.15, to=-0.5, by=-0.05)
		Dup <- seq(from=0.15, to=0.5, by=0.05)
		ChrMean <- seq(from=-0.4, to=0.2, by=0.025)
		#				                                                                 0  .25 .05 .075 .1  .125 .15  .175 0.2	
		ChrMeanProb <- c(0.01, 0.01,0.01,0.01,0.01, 0.01,0.01,0.01,0.02,0.02,0.1,0.15,0.2,0.1,0.4,0.6,0.8,0.6,0.4,0.1,0.02,0.01,0.01,0.01,0.01)

		ChrSD <- seq(from=0.1, to=0.5, by=0.025)
		#	       0.1 .125  .15 .175  0.2  .225 .25  .275  .3 .325 .35 .375 .4 .425 .45 .475 .5 	
		ChrSDProb <- c(0.2, 0.4, 0.6, 0.8, 0.6,  0.4, 0.3, 0.2, 0.3, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01)
	
		# BAF Del prob
		BAF_Del <- BAF_Basic + 0.02
		BAF_Del[1:2] <- BAF_Del[1:2] + 0.7
		BAF_Del[100:101] <- BAF_Del[100:101] + 0.9
		
		# BAF Dup prob
		BAF_Dup <- BAF_Basic + 0.02
		BAF_Dup[1:2] <- BAF_Dup[1:2] + 0.4
		BAF_Dup[100:101] <- BAF_Dup[100:101] + 0.6
		BAF_Dup[25:35] <- BAF_Dup[25:35] + 0.15 # 0.25 0.30 0.35, 
		BAF_Dup[70:80] <- BAF_Dup[70:80] + 0.15 # 0.7 0.75 0.8
		
	}
	
	tmp <- list(Del=Del, Dup=Dup, ChrMean=ChrMean, ChrSD=ChrSD, ChrSDProb=ChrSDProb, TelomereNoiseSize=TelomereNoiseSize, TelomereNoiseEffect=TelomereNoiseEffect, BAF_Normal=BAF_Normal, BAF_Del=BAF_Del, BAF_Dup=BAF_Dup, BadSNPs=BadSNPs, BadSNPIntensity=BadSNPIntensity, BadSNPIntensityProb=BadSNPIntensityProb, BAFs=BAFs, ChrMeanProb=ChrMeanProb)
	return(tmp)
}
