##' FindCNV: Find CNVs using Log.R.Ratio data by cpt.meanvar method. 
##'
##' Input: sample path, CNV in data frame format and min number of SNPs
##' @title FindCNV
##' @return A data frame with CNVs, with start and stop position, chr and etc. for each sample. 
##' @author Marcelo Bertalan
##' @export

FindCNV.V4 <- function(ID, MinNumSNPs, CNV)
{
	
	library("changepoint")

	tmp <- sapply(unique(CNV$Chr), function(X) # X is chr. Loop over Chr.
	{
		CHR <- as.character(X)
		subCNV <- subset(CNV, Chr %in% CHR)
		subCNV <- subCNV[with(subCNV, order(Position)),]
		
		# Using changepoint package	
		CPT.Res <- cpt.mean(subCNV$Log.R.Ratio, method='PELT', penalty="AIC") 
		indx <- cpts(CPT.Res)
		indx <- c(1, indx, length(subCNV$Log.R.Ratio))


		# Mergin CNVs #
		#Info <- sapply(1:(length(indx)-1), function(X){ rbind(indx[X],indx[(X+1)]) })
		#Info <- t(Info)
		#Df <- as.data.frame(Info)
		#CNVMean <- apply(Df, 1, function(X){ mean(CPT.Res@data.set[X[1]:X[2]]) } )
		#Df$CNVMean <- CNVMean
		#names(Df) <- c("Start", "Stop", "CNVMean")
		#tmp <- apply(Df, 1, function(X){ rep(X[3], (as.numeric(X[2])-as.numeric(X[1]))) })
		#tmp <- unlist(tmp)
		#tmp[tmp < -0.15] <- -0.5
		#tmp[tmp > 0.15] <- 0.5
		#tmp <- as.vector(tmp)
		#if(length(tmp) > 10)
		#{
		#	CPT.Res.tmp <- cpt.mean(tmp, method='PELT', penalty="AIC")
		#	indx <- cpts(CPT.Res.tmp)
		#	indx <- c(1, indx, length(tmp))
		#}
		###############

		DF <- DefineStartAndStop(indx, subCNV, MinNumSNPs, CHR, ID, CPT.Res)
		Df <- FixCNVPosition(DF, subCNV, MinNumSNPs, ID) # Fix the position of the CNV by distance of the SNPs
		return(Df)
	})

	df <- MatrixOrList2df(tmp)
	if(nrow(df) > 0)
	{
		df <- subset(df, NumSNPs > MinNumSNPs)
		df <- df[,!colnames(df) %in% ".id"]
		return(df)
	}
}

