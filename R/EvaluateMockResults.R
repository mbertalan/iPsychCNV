##' EvaluateMockResults: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title EvaluateMockResults
##' @param MockCNV: Unknown. 
##' @param Df: Data frame with predicted CNVs for each sample? Unknown.
##' @param Cores: Number of cores used, default = 1.
##' @param MinOverlap, Minimum overlap to accept a CNV, default = 80%.
##' @param MaxOverlap, maximum overlap to accept a CNV, default = 120%.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

EvaluateMockResults <- function(MockCNVs, df, Cores=1, MinOverlap=80, MaxOverlap=120)
{
	library(pROC)
	library(parallel)
	#if(length(MockCNVs$CNVID) == 0){ MockCNVs$CNVID <- 1:nrow(MockCNVs) } # if no CNVID is provide
	#if(length(df$CNVID) == 0){ df$CNVID <- 1:nrow(df) } # if no CNVID is provide
	#if(length(MockCNVs$ID) == 0){ MockCNVs$ID <- "Sample" } # if no CNVID is provide
	#if(length(df$ID) == 0){ df$ID <- "Sample" } # if no CNVID is provide
	
	if(length(df$probs) == 0){ df$probs <- 1 }
	
	df <- df[, c("CN", "Start", "Stop", "CNVID", "ID", "Length", "Chr", "NumSNPs", "probs")]
	df$Start <- as.numeric(df$Start)
	df$Stop <- as.numeric(df$Stop)
	df$Chr <- as.numeric(df$Chr)
	df$CN <- as.numeric(df$CN)
	df$Length <- as.numeric(df$Length)
	df$NumSNPs <- as.numeric(df$NumSNPs)
	df$probs <- as.numeric(df$probs)
	
	Eval <- mclapply(1:nrow(MockCNVs), mc.cores=Cores, mc.preschedule = FALSE, function(i)
	{
		X <- MockCNVs[i,]
		StartM <- as.numeric(X["Start"])
		StopM <- as.numeric(X["Stop"])
		IDM <- as.character(X["ID"])
		ChrM <- as.numeric(X["Chr"])
		LengthM <- as.numeric(X["Length"])
		NumSNPsM <- as.numeric(X["NumSNPs"])
		CNVID <- as.numeric(X["CNVID"])
		CNM <- as.numeric(X["CN"])

		if(CNM == 2){ CNV.Present <- 0 }else{ CNV.Present <- 1 }
		# From df predition
 		res <- subset(df, Chr == ChrM & Start <= StopM & Stop >= StartM & ID %in% IDM)
 	
		NumCNVs <- nrow(res)
		if(NumCNVs > 0) # Selecting only one CNV
		{
			DifLength <- abs(res$Length - LengthM)
			names(DifLength) <- res$Length
			res <- subset(res, Length == as.numeric(names(sort(DifLength))[1]))
			#res <- subset(res, Length == max(res$Length))	
			if(nrow(res) > 1){ res <- res[1,]}		
			CNVID2 <- res$CNVID
			CN2 <- as.numeric(res$CN)
			OneProb <- as.numeric(res$probs)
		}

		if(NumCNVs == 0)
		{
			CNV.Predicted <- 0
			PredictedByOverlap <- 0
			OverlapLenghM <- 0
			OverlapSNP <-0
			CNVID2 <- NA
			CN2 <- 2
			OneProb <- 0
			
		}
		else
		{	
			StartO <- max(c(as.numeric(res$Start), as.numeric(StartM)))[1]
			StopO <- min(c(as.numeric(res$Stop), as.numeric(StopM)))[1]
			if(StartO == StartM & StopO == StopM)
			{
				StartO <- max(res$Start)[1]
				StopO <- min(res$Stop)[1]
			}
			LengthO <- StopO - StartO
			OverlapLenghM <- ((LengthO/LengthM)*100)
			OverlapSNP <- ((res$NumSNPs/NumSNPsM)*100)

			if(CNM == 2)
			{
				res2 <- subset(df, Chr == ChrM &  StartM > Stop & StopM > Start & ID %in% IDM) # If is a non-CNV region, any overlap will count.				if(nrow(res2) > 0)
				{
					CNV.Predicted <- OneProb
					PredictedByOverlap <- 0		# It doesn't matter the overlap the prediction is in a non-CNV region.
				}
				else
				{
					CNV.Predicted <- 0
					PredictedByOverlap <- 0
				}
			}
			else
			{
				if(CNM == CN2 & CNM != 2)
				{	
					CNV.Predicted <- OneProb
					if(OverlapLenghM > MinOverlap & OverlapLenghM < 120)
					{ 
						PredictedByOverlap <- OverlapLenghM
					}
					else
					{
						PredictedByOverlap <- 0 	# CNV predicted is far from what it should be
					}
				}
				else						# There is a prediction but is not correct
				{ 
					CNV.Predicted <- OneProb
					PredictedByOverlap <- 0
				} 
			}
		}
		cat("ChrM:",ChrM, "CNM", CNM, "CN2", CN2, "CNV.Predicted", CNV.Predicted, "OverlapLenghM", OverlapLenghM, "OverlapSNP", OverlapSNP, "CNVID",CNVID, "CNVID2",CNVID2, "CN2", CN2, "NumCNVs", NumCNVs, "PredictedByOverlap", PredictedByOverlap,"\r") 
		
		df2 <- data.frame(CNV.Present=CNV.Present, CNV.Predicted=CNV.Predicted, Overlap.Length=OverlapLenghM, Overlap.SNP=OverlapSNP, CNVID.Mock=CNVID, CNVID.Pred=CNVID2, CN.Pred=CN2, NumCNVs = NumCNVs, PredictedByOverlap=PredictedByOverlap, stringsAsFactors=FALSE)
		save(df2, file="df2.RData")
		return(df2)
	})
	#save(Eval, file="Eval.RData")
	Eval <- MatrixOrList2df(Eval)
	
	if(nrow(Eval) == nrow(MockCNVs))
	{
		df3 <- cbind(Eval, MockCNVs)
	}
	else
	{
		df3 <- Eval
		stop("More samples than CNV or vice versa. Evaluation to do match MockCNVs.\n")
	}
	return(df3)
}
