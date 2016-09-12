##' GetCNVFreq: Get CNV freq. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetCNVFreq
##' @param df Data frame with CNVs, must have status1 as phenotype variable and ID as sample ID variable.
##' @param Pheno Data frame with phenotypes must have status1 as phenotype variable and ID as sample ID variable.
##' @param Cut Cutoff value for CNV mean. Only negative for now. Default Cuf=-0.3.
##' @param roi Regions of interest, object with information where should be the CNV you are looking for. 
##' @return return data in a data frame.
##' @author Marcelo Bertalan.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

GetCNVFreq <- function(df, Pheno, roi, OverlapMin=0.9, OverlapMax=1.1)
{
	Source <- unique(df$Source)
	TotalPheno <- table(Pheno$status1)
	tmp <- SelectSamplesFromROI(DF=df, roi=roi, OverlapMin=OverlapMin, OverlapMax=OverlapMax)
	#tmp <- subset(tmp, CNVmean < Cut)
	cat(roi$locus, "\n")
	if(nrow(tmp) > 1)
	{
		tmp <- tmp[!duplicated(tmp[,"ID"]),]
		TotalCNVs <- table(tmp$status1)
	
		indx <- intersect(names(TotalCNVs), names(TotalPheno))
		Perc <- (TotalCNVs[indx]/TotalPheno[indx])*100

		df <- data.frame(status1=indx, CNV=as.vector(TotalCNVs[indx]), Samples=as.vector(TotalPheno[indx]), Freq=as.vector(Perc), Source=Source, Cut=Cut, Chr=roi$Chr, Start=roi$Start, Stop=roi$Stop, locus=roi$locus, stringsAsFactors=FALSE)
		return(df)
	}
}
