##' ReadCNV: read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ReadCNV
##' @return return data in data frame
##' @author Marcelo Bertalan
##' @export

ReadCNV <- function(RawFile="Test.txt", skip=0)
{
	library(data.table)
	CNV <- fread(RawFile, head=T, sep="\t", skip=skip)
	CNV <- as.data.frame(CNV)
	colnames(CNV) <- gsub(" ", ".", colnames(CNV))
	CNV <- CNV[,c("Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq")] # SNP.Name
	CNV <- subset(CNV, !Chr %in% c("MT", "X", "Y", "XY", "0"))
	CNV <- subset(CNV, !is.na(CNV$B.Allele.Freq))
	CNV <- subset(CNV, !is.na(CNV$Log.R.Ratio))
	BadSNPs <- c("psy_rs174693","exm1585879","rs2239393","kgp11396666","kgp15052777","exm1585985","exm1999632","exm1586022","exm2276749","exm1586057","rs887204","exm1586114","exm1586120","rs756653","rs2073746")
	#CNV <- subset(CNV, !SNP.Name %in% BadSNPs) 
	CNV$LRR <- CNV$Log.R.Ratio # CNV$LRR is the
	return(CNV)
}
