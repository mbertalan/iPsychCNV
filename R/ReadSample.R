##' ReadSample: read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ReadSample
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

ReadSample <- function(RawFile="Test.txt", skip=0, LCR=FALSE, PFB=NULL, chr=NA)
{
	suppressPackageStartupMessages(library(data.table))
	CNV <- fread(RawFile, head=T, sep="\t", skip=skip, showProgress=FALSE, verbose=FALSE)
	CNV <- as.data.frame(CNV)
	colnames(CNV) <- gsub(" ", ".", colnames(CNV))
	colnames(CNV)[colnames(CNV) %in% "Name"] <- "SNP.Name"
	colnames(CNV)[colnames(CNV) %in% "Chromosome"] <- "Chr"
	colnames(CNV)[colnames(CNV) %in% "Allele1.-.Top"] <- "Allele1"
	colnames(CNV)[colnames(CNV) %in% "Allele2.-.Top"] <- "Allele2"
	colnames(CNV)[colnames(CNV) %in% "BAF"] <- "B.Allele.Freq"
	colnames(CNV)[colnames(CNV) %in% "LRR"] <- "Log.R.Ratio"
	#CNV <- CNV[,c("SNP.Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "Allele1", "Allele2")] # SNP.Name
	
	# removing chr from chromosome name (deCODE)
	CNV$Chr <- gsub("chr", "", CNV$Chr)
	
	# Genotype together (deCODE)
	if(!is.null(CNV$Genotype))
	{
		Allele <- sapply(CNV$Genotype, function(X){ data.frame(Allele1= unlist(strsplit(X, ""))[1], Allele2=unlist(strsplit(X, ""))[2], stringsAsFactors=F) })
		CNV$Allele1 <- Allele$Allele1
		CNV$Allele2 <- Allele$Allele2
	}
	
	# PFB
	if(is.null(PFB)){ CNV$PFB <- rep(0.5, nrow(CNV)) }else{ CNV$PFB <- PFB }
	
	# Subseting 
	CNV <- subset(CNV, !Chr %in% c("XY", "0")) # "MT", "X", "Y",
	CNV <- subset(CNV, !is.na(CNV$B.Allele.Freq))
	CNV <- subset(CNV, !is.na(CNV$Log.R.Ratio))
	
	# chr specific. Example chr="22"
	if(!is.na(chr)){ CNV <- subset(CNV, Chr %in% chr) }
	
	if(LCR == TRUE)
	{
		CNV <- subset(CNV, !SNP.Name %in% LCR.SNPs) 
	}
	CNV$LRR <- CNV$Log.R.Ratio # CNV$LRR is the
	return(CNV)
}
