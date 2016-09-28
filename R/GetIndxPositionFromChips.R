##' GetIndxPositionFromChips: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetIndxPositionFromChips
##' @param CNVs: Data frame with information on CNVs. Unknown.
##' @param Sample: Unknown.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

GetIndxPositionFromChips <- function(CNVs, Sample)
{
	Res2 <- sapply(unique(CNVs$Chr), function(y)
	{
		CHR2 <- y
		CHR2 <- gsub(" ", "", CHR2)
		df2 <- subset(CNVs, Chr %in% CHR2)	
		subSample <- subset(Sample, Chr %in% CHR2) 			# Subset by Chr
		subSample <- subSample[with(subSample, order(Position)),]	
	
		Res <- sapply(1:nrow(df2), function(i)
		{
			X <- df2[i,]
			StartM <- as.numeric(X["Start"]) # CNV Start
			StopM <- as.numeric(X["Stop"])
			# Order
			StartIndx <- which.min(abs(subSample$Position - StartM))
			StopIndx <- which.min(abs(subSample$Position - StopM))
			NumSNPs <- StopIndx - StartIndx
			
			StartM <- subSample$Position[StartIndx]
			StopM <- subSample$Position[StopIndx]
			
			SNP.Start <- subSample$SNP.Name[StartIndx]
			SNP.Stop <- subSample$SNP.Name[StopIndx]
			
			df <- data.frame(Start=StartM, Stop=StopM, Chr=CHR2, StartIndx=StartIndx, StopIndx=StopIndx, NumSNPs=NumSNPs, SNP.Start=SNP.Start, SNP.Stop=SNP.Stop, stringsAsFactors=F)
			#df <- data.frame(StartIndx=StartIndx, StopIndx=StopIndx, NumSNPs=NumSNPs)
			df3 <- cbind(df, df2[i,])
			return(df3)
		})
		df <- MatrixOrList2df(tmp=Res)
		return(df)
	})
	df <- MatrixOrList2df(tmp=Res2)
	return(df)
}
