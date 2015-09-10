##' MultipleMockData: Runs mock data multiple times. 
##'
##' Specifically designed to handle noisy data from amplified DNA on dried blood spot (DBS) cards. The function is a pipeline using many subfunctions.
##' @title MultipleMockData
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

MultipleMockData <- function(NSamples=10, NLoops=10, Cores=30)
{
	Res <- sapply(1:NLoops, function(Loops)
	{
		# Creating Mock data
		MockDataCNVs <- MockData(N=NSamples, Type="PKU", Cores=Cores)
	
		# iPsychCNV prediction
		iPsych.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=20, Cores=Cores, Pattern="^MockSample", MinLength=10, Skip=0, LCR=FALSE, Quantile=FALSE)

		# PennCNV
		PennCNV.Pred <- RunPennCNV(PathRawData=".", Pattern="^MockSample_*", Cores=Cores, Skip=0, Normalization=FALSE, PFB=NA, HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/" )

		# Evaluating methods
		iPsychCNV.Eval <- EvaluateMockResults(MockDataCNVs, iPsych.Pred)
		PennCNV.Eval <- EvaluateMockResults(MockDataCNVs, PennCNV.Pred)

		# When prediction fail.
		## PennCNV
		tmp <- subset(PennCNV.Eval, CN != 2)
		tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })

		df <- data.frame(True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="PennCNV", Loop=Loops, stringsAsFactors=F)

		## iPsychCNV
		tmp <- subset(iPsychCNV.Eval, CN != 2)
		tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })

		df2 <- data.frame(True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="iPsychCNV", Loop=Loops, stringsAsFactors=F)
		df3 <- rbind(df, df2)
		system("rm -f MockSample_*")
		return(df3)
	})
	Res <- MatrixOrList2df(Res)
	Res$Values <- as.numeric(Res$Values)
	return(Res)
}
# Res <- MultipleMockData(NSamples=30, NLoops=30)
# ggplot(Res, aes(x=as.factor(Values), y=True.positive)) + geom_boxplot(aes(fill=as.factor(Source)), outlier.shape=NA) +  facet_wrap(~ By, scales = "free") + scale_fill_brewer(palette="Set1")
