##' MultipleMockData: Runs mock data multiple times. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MultipleMockData
##' @param NSamples: Unknown, default = 10.
##' @param NLoops: Unknown, default = 10.
##' @param Cores: Number of cores used, default = 30.
##' @param HMM: Unknown, default = Unknown.
##' @param Path2PennCNV: The path to pennCNV.
##' @return Return data in a data frame.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

MultipleMockData <- function(NSamples=100, NLoops=10, Cores=28, HMM="/media/NeoScreen/NeSc_home/share/Programs/penncnv/lib/hhall.hmm", Path2PennCNV="/media/NeoScreen/NeSc_home/share/Programs/penncnv/")
{
	library(pROC)
	Res <- sapply(1:NLoops, function(Loops)
	{
		# Creating Mock data
		cat("Creating mock data: Loop ", Loops, "\n")
		MockDataCNVs <- MockData(N=NSamples, Type="PKU", Cores=Cores)
		MockHot <- subset(MockDataCNVs, ID %in% "MockSample_1.tab" & CN != 2)
		
		# iPsychCNV prediction
		cat("Running iPsychCNV: Loop ", Loops, "\n" )
		iPsych.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=20, Cores=Cores, Pattern="^MockSample", MinLength=10, Skip=0, LCR=FALSE, Quantile=FALSE)
		iPsych.Pred <- subset(iPsych.Pred, CN != 2)
		
		# iPsychCNV + Hotspots + ReScanCNVs
		#CNVs.Hotspots <- HotspotsCNV(df=iPsych.Pred, Freq=3, OverlapCutoff=0.8, Cores=Cores)
		#iPsych.rescan <- ReScanCNVs(CNVs=MockHot, Cores=Cores, Pattern="^MockSample_*", Skip=0, hg="hg19", PathRawData=".", IndxPos=TRUE)
	
		# PennCNV
		cat("Running PennCNV: Loop ", Loops, "\n")
		PennCNV.Pred <- RunPennCNV(PathRawData=".", Pattern="^MockSample.*", Cores=Cores, Skip=0, Normalization=FALSE, PFB="NO", HMM=HMM, Path2PennCNV=Path2PennCNV)

		# PennCNV Filter
		#PennCNV.filter <- FilterFromCNVs(CNVs=PennCNV.Pred, PathRawData=".", MinNumSNPs=10, Source="PennCNV.Filter", Skip=0, Cores=Cores)	

		# Filter + Hotspots + ReScanCNVs
		#tmp <- subset(PennCNV.filter, CN != 2)
		#CNVs.Hotspots <- HotspotsCNV(df=tmp[,1:20], Freq=1, OverlapCutoff=0.8, Cores=Cores)
		#PennCNV.rescan <- ReScanCNVs(CNVs=CNVs.Hotspots, Cores=Cores, Pattern="^MockSample_*", Skip=0, hg="hg19", PathRawData=".")
		#PennCNV.rescan$ID <- PennCNV.rescan$SampleID

		# Evaluating methods
		cat("iPsychCNV evaluting methods: Loop ", Loops, "\n")
		iPsychCNV.Eval <- EvaluateMockResults(MockDataCNVs, iPsych.Pred, Cores=Cores)
		cat("PennCNV Evaluting methods: Loop ", Loops, "\n")
		PennCNV.Eval <- EvaluateMockResults(MockDataCNVs, PennCNV.Pred, Cores=Cores)
		#cat("ReiPsychCNV Evaluting methods: Loop ", Loops, "\n")
		#ReiPsych.Eval <- EvaluateMockResults(MockDataCNVs, iPsych.rescan, Cores=Cores)
		cat("Done evaluting methods: Loop ", Loops, "\n")
		#Filter.Eval <- EvaluateMockResults(MockDataCNVs, PennCNV.filter, Cores=Cores)
		#Rescan.Eval <- EvaluateMockResults(MockDataCNVs, PennCNV.rescan, Cores=Cores)
		
		# ROC 
		iPsychCNV.AUC <- roc(iPsychCNV.Eval$CNV.Predicted, iPsychCNV.Eval$CNV.Present)
		#ReiPsych.AUC <- roc(ReiPsych.Eval$CNV.Predicted, ReiPsych.Eval$CNV.Present)
		PennCNV.AUC <- roc(PennCNV.Eval$CNV.Predicted, PennCNV.Eval$CNV.Present)
		#Filter.AUC <- roc(Filter.Eval$CNV.Predicted, Filter.Eval$CNV.Present)$auc[1]
		#Rescan.AUC <- roc(Rescan.Eval$CNV.Predicted, Rescan.Eval$CNV.Present)$auc[1]

		# Coords
		iPsychCNV.AUC.coords <- coords(iPsychCNV.AUC, "best", ret=c("specificity", "sensitivity", "accuracy", "tn", "fn", "fp", "npv", "ppv"))
		#ReiPsych.AUC.coords <- coords(ReiPsych.AUC, "best", ret=c("specificity", "sensitivity", "accuracy", "tn", "fn", "fp", "npv", "ppv"))
		PennCNV.AUC.coords <- coords(PennCNV.AUC, "best", ret=c("specificity", "sensitivity", "accuracy", "tn", "fn", "fp", "npv", "ppv"))
		
		iPsychCNV.AUC.coords <- as.data.frame(t(iPsychCNV.AUC.coords))
		#ReiPsych.AUC.coords <- as.data.frame(t(ReiPsych.AUC.coords))
		PennCNV.AUC.coords <- as.data.frame(t(PennCNV.AUC.coords))
		
		iPsychCNV.AUC.coords$auc <- iPsychCNV.AUC$auc[1]
		#ReiPsych.AUC.coords$auc <- ReiPsych.AUC$auc[1]
		PennCNV.AUC.coords$auc <-  PennCNV.AUC$auc[1]
		
		iPsychCNV.AUC.coords$source <- "iPsychCNV"
		#ReiPsych.AUC.coords$source <- "Re-ScanCNVs"
		PennCNV.AUC.coords$source <- "PennCNV"
		
		# When prediction fail.
		## PennCNV
		tmp <- subset(PennCNV.Eval, CN != 2)
		tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })
		tmp6 <- tapply(tmp$True.Positive, as.factor(tmp$Heterozygosity), function(X){ sum(X)/length(X) })
		df <- data.frame(AUC=PennCNV.AUC$auc[1], True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5), as.numeric(tmp6)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5)),  rep("Heterozygosity", length(tmp6))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5), names(tmp6)), Source="PennCNV", Loop=Loops, stringsAsFactors=F)

		## iPsychCNV
		tmp <- subset(iPsychCNV.Eval, CN != 2)
		tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })
		tmp6 <- tapply(tmp$True.Positive, as.factor(tmp$Heterozygosity), function(X){ sum(X)/length(X) })
		#df2 <- data.frame(AUC=iPsychCNV.AUC,True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="iPsychCNV", Loop=Loops, stringsAsFactors=F)
		df2 <- data.frame(AUC=iPsychCNV.AUC$auc[1], True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5), as.numeric(tmp6)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5)),  rep("Heterozygosity", length(tmp6))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5), names(tmp6)), Source="iPsychCNV", Loop=Loops, stringsAsFactors=F)

		## iPsychCNV + ReScan
		#tmp <- subset(ReiPsych.Eval, CN != 2)
		#tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		#tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		#tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		#tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		#tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })
		#tmp6 <- tapply(tmp$True.Positive, as.factor(tmp$Heterozygosity), function(X){ sum(X)/length(X) })
		#df3 <- data.frame(AUC=ReiPsych.AUC,True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="iPsychCNV + ReScan", Loop=Loops, stringsAsFactors=F)
		#df3 <- data.frame(AUC=ReiPsych.AUC, True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5), as.numeric(tmp6)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5)),  rep("Heterozygosity", length(tmp6))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5), names(tmp6)), Source="iPsychCNV + ReScan", Loop=Loops, stringsAsFactors=F)

		## PennCNV + Filter
		#tmp <- subset(Filter.Eval, CN != 2)
		#tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		#tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		#tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		#tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		#tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })
		#df3 <- data.frame(AUC=Filter.AUC, True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="Filter", Loop=Loops, stringsAsFactors=F)
		
		## PennCNV + Filter + Hotspots + ReScan
		#tmp <- subset(Rescan.Eval, CN != 2)
		#tmp$True.Positive <- tmp$CNV.Present == tmp$CNV.Predicted
		#tmp2 <- tapply(tmp$True.Positive, as.factor(tmp$sd), function(X){ sum(X)/length(X) })
		#tmp3 <- tapply(tmp$True.Positive, as.factor(tmp$CNVmean), function(X){ sum(X)/length(X) })
		#tmp4 <- tapply(tmp$True.Positive, as.factor(tmp$NumSNPs), function(X){ sum(X)/length(X) })
		#tmp5 <- tapply(tmp$True.Positive, as.factor(tmp$CN), function(X){ sum(X)/length(X) })
		#df4 <- data.frame(AUC=Rescan.AUC, True.positive=c(as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5)), By=c(rep("sd", length(tmp2)), rep("CNV.mean", length(tmp3)), rep("NumSNPs", length(tmp4)), rep("CN", length(tmp5))), Values=c(names(tmp2), names(tmp3), names(tmp4), names(tmp5)), Source="ReScan", Loop=Loops, stringsAsFactors=F)
		
		df5 <- rbind(df, df2)
		#df5 <- rbind(df, df2, df3)
		#df5 <- rbind(df, df2, df3, df4)
		system("rm -f MockSample_*")
		return(df5)
	})
	Res <- MatrixOrList2df(Res)
	Res$Values <- as.numeric(Res$Values)
	return(Res)
}
# Res <- MultipleMockData(NSamples=30, NLoops=30)
# ggplot(Res, aes(x=as.factor(Values), y=True.positive)) + geom_boxplot(aes(fill=as.factor(Source)), outlier.shape=NA) +  facet_wrap(~ By, scales = "free") + scale_fill_brewer(palette="Set1")
