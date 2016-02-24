##' RunGada: Unknown.
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title RunGada
##' @param MockData: The mock data set, default = tmp.
##' @param MinSegLen: Unknown, default = 50.
##' @param Tgada: Unknown, default = 4.5.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

RunGada <- function(MockData=tmp, MinSegLen=50, Tgada=4.5)
{
	library(gada)
	# Creating gen.info for Gada
	gen.info <- MockData[,1:3]
	colnames(gen.info) <- c("Name", "chr", "pos")
	
	# Gada can not run just one chromosome.
	NotUsed <- gen.info
	NotUsed$chr <- 2
	gen.info <- rbind(gen.info, NotUsed)
	
	dataMock<-setupGADAgeneral(c(MockData$Log.R.Ratio,MockData$Log.R.Ratio) ,gen.info=gen.info)
	step1<-SBL(dataMock,  estim.sigma2=TRUE)
	step2<-BackwardElimination(step1,T=Tgada,MinSegLen=MinSegLen)
	tmp <- summary(step2)
	tmp <- tmp[,1:6]
	colnames(tmp) <- c("Start","Stop", "NumSNPs", "CNVMean", "Chr", "CN")
	tmp <- subset(tmp, Chr %in% "1")
	tmp$CN[tmp$CN == 0] <- 2
	tmp$CN[tmp$CN == 1] <- 3
	tmp$CN[tmp$CN == -1] <- 1
	tmp$ID <- rep("Mock", nrow(tmp)) 
	tmp$Source <- "Gada"
	tmp$Length <- tmp$Stop - tmp$Start
	tmp$CNVID <- 1:nrow(tmp)
	PredictedCNV <- tmp
	return(PredictedCNV)
}
