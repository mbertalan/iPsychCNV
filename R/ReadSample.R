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
	Sample <- fread(RawFile, head=T, sep="\t", skip=skip, verbose=FALSE)
	Sample <- as.data.frame(Sample)
	colnames(Sample) <- gsub(" ", ".", colnames(Sample))
	colnames(Sample)[colnames(Sample) %in% "Name"] <- "SNP.Name"
	colnames(Sample)[colnames(Sample) %in% "Chromosome"] <- "Chr"
	colnames(Sample)[colnames(Sample) %in% "Allele1.-.Top"] <- "Allele1"
	colnames(Sample)[colnames(Sample) %in% "Allele2.-.Top"] <- "Allele2"
	colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
	colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
	#CNV <- CNV[,c("SNP.Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "Allele1", "Allele2")] # SNP.Name
	
	# removing chr from chromosome name (deCODE)
	Sample$Chr <- gsub("chr", "", Sample$Chr)
	
	# Genotype together (deCODE)
	#if(!is.null(Sample$Genotype))
	#{
	#	Allele <- sapply(Sample$Genotype, function(X){ data.frame(Allele1= unlist(strsplit(X, ""))[1], Allele2=unlist(strsplit(X, ""))[2], stringsAsFactors=F) })
	#	Allele <- MatrixOrList2df(Allele)
	#	Sample$Allele1 <- Allele$Allele1
	#	Sample$Allele2 <- Allele$Allele2
	#}
	
	# PFB
	if(is.null(PFB)){ Sample$PFB <- rep(0.5, nrow(Sample)) }else{ Sample$PFB <- PFB }
	
	# Subseting 
	Sample <- subset(Sample, !Chr %in% c("XY", "0")) # "MT", "X", "Y",
	Sample <- subset(Sample, !is.na(Sample$B.Allele.Freq))
	Sample <- subset(Sample, !is.na(Sample$Log.R.Ratio))
	
	# chr specific. Example chr="22"
	if(!is.na(chr)){ Sample <- subset(Sample, Chr %in% chr) }
	
	if(LCR == TRUE)
	{
		Sample <- subset(Sample, !SNP.Name %in% LCR.SNPs) 
	}
	
	Sample$LRR <- Sample$Log.R.Ratio # CNV$LRR is the
	return(Sample)
}
