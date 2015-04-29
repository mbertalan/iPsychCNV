RunGada <- function(MockData=tmp)
{
	library(gada)
	# Creating gen.info for Gada
	gen.info <- MockData[,1:3]
	colnames(gen.info) <- c("Name", "chr", "pos")
	Start <- round(nrow(MockData)/2) # Gada can not run just one chromosome.
	gen.info$chr[Start:nrow(MockData)] <- "2"
	dataMock<-setupGADAgeneral(MockData$Log.R.Ratio,gen.info=gen.info)
	step1<-SBL(dataMock,  estim.sigma2=TRUE)
	step2<-BackwardElimination(step1,T=4.5,MinSegLen=3)
	tmp <- summary(step2)
	tmp <- tmp[,1:6]
	colnames(tmp) <- c("Start","Stop", "NumSNPs", "CNVMean", "Chr", "CN")
	tmp$CN[tmp$CN == 0] <- 2
	tmp$CN[tmp$CN == 1] <- 3
	tmp$CN[tmp$CN == -1] <- 1
	tmp$Chr <- rep("1", nrow(tmp))
	tmp$ID <- rep("Mock", nrow(tmp)) 
	PredictedCNV <- tmp
	return(PredictedCNV)
}
