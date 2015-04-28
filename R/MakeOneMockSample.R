##' MakeOneMockSample: Run a single mock sample   
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MakeOneMockSample
##' @return Classification for LRR.
##' @author Marcelo Bertalan
##' @export

MakeOneMockSample <- function(Noise=0.5, CNVMean=0.2)
{
	library(RColorBrewer)
	library(ggplot2)
	library(ggbio)

	CNVs_Mean <- rep(CNVMean, 24)
	names(CNVs_Mean) <- c(3, 1, 1, 3, 0, 3, 4, 3, 3, 3, 4, 3, 3, 1, 4, 1, 1, 0, 3, 1, 1, 1, 0, 1)
	CNVs_Mean[names(CNVs_Mean) %in% "1"] <- CNVs_Mean[names(CNVs_Mean) %in% "1"] * -1
	CNVs_Mean[names(CNVs_Mean) %in% "0"] <- CNVs_Mean[names(CNVs_Mean) %in% "0"] * -1
	
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

	# BAF CN=0
	BAF_CN0 <- BAF_Basic
	
	# BAF Dup prob
	BAF_Dup <-  BAF_Basic
	BAF_Dup[1:2] <- BAF_Dup[1:2] + 1
	BAF_Dup[100:101] <- BAF_Dup[100:101] + 1
	BAF_Dup[30:35] <- BAF_Dup[30:35] + 0.05 
	BAF_Dup[32:33] <- BAF_Dup[32:33] + 0.1 
	BAF_Dup[65:70] <- BAF_Dup[65:70] + 0.05
	BAF_Dup[67:68] <- BAF_Dup[67:68] + 0.1 	

	# BAF CN=4
	BAF_CN4 <- BAF_Dup
	BAF_CN4[c(47:53)] <- BAF_CN4[c(47:53)] + 0.05
	BAF_CN4[c(50:51)] <- BAF_CN4[c(50:51)] + 0.1

	BAF <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=50000)
	LRR <- rnorm(50000, mean=0, sd=0.1)

	sapply(1:24, function(i)
	{
		Start <- i * 2000
		Stop <- Start+500
		CNVmean <- CNVs_Mean[i] # sample(c(CNVMean, (CNVMean*-1)), 1)
		LRR[Start:Stop] <<- LRR[Start:Stop] + CNVmean
		if(names(CNVmean) %in% "3")
		{
			BAFCNV <- sample(BAFs, prob=BAF_Dup, replace=TRUE, size=501)
			BAF[Start:Stop] <<- BAFCNV
		}
		else if(names(CNVmean) %in% "1")
		{	
			BAFCNV <- sample(BAFs, prob=BAF_Del, replace=TRUE, size=501)
			BAF[Start:Stop] <<- BAFCNV
		}
		else if(names(CNVmean) %in% "0")
		{
			BAFCNV <- sample(BAFs, prob=BAF_CN0, replace=TRUE, size=501)
			BAF[Start:Stop] <<- BAFCNV
		}
		else if(names(CNVmean) %in% "4")
		{
			BAFCNV <- sample(BAFs, prob=BAF_CN4, replace=TRUE, size=501)
			BAF[Start:Stop] <<- BAFCNV
		}
	})		

	CNV.chr1 <- subset(CNV, Chr %in% "1")
	SNP.Name <- CNV.chr1$SNP.Name[1:50000]
	Chr <- CNV.chr1$Chr[1:50000]
	Position <- CNV.chr1$Position[1:50000]

	if(Noise == 1) # Random BAF. CN = 4
	{
		BAF <- sample(BAFs, prob=BAF_Basic, replace=TRUE, size=50000)
	}
	if(Noise == 2) # No heterozygosity. CN = 1
	{
		BAF <- sample(BAFs, prob=BAF_Del, replace=TRUE, size=50000)
	}
	if(Noise == 3) #Only heterozygosity. CN = 2
	{
		BAF <- sample(BAFs, prob=BAF_Normal, replace=TRUE, size=50000)
	}


	df <- data.frame(SNP.Name=SNP.Name, Chr=Chr, Position=Position, Log.R.Ratio=LRR, B.Allele.Freq=BAF, stringsAsFactors=F)
	df$CNVmean <- CNVMean
	df$CNVmean[df$CN == 1] <- df$CNVmean[df$CN == 1] * -1
	write.table(df, sep="\t", quote=FALSE, row.names=FALSE, file="MockSample_1.tab") 
	return(df)
}
