##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays.
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title iPsychCNV
##' @param Files: Object with full path to all Illumina files.
##' @param Size: Number of files.
##' @param SNPList: Getting Chr. and Position from another source than the RawFile - input should be the full path of the SNPList with columns: Name, Chr, amd Position. Any positions from the RawFile will be erased. A PFB-column is also allowed but will be overwritten by the PFB-parameter or exchanged with 0.5
##' @param Cores: Number of cpus.
##' @param Skip: Number of lines to be skipped.
##' @param chr: specify which chromosomes should be analyzed. 
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export

GetHeterozygosity <- function(Files=Files, Size=1000, Cores=28, Skip=10, SNPList=NULL, chr=1:22)
{
	suppressPackageStartupMessages(library(parallel))
	Indx <- sample(x=length(Files), size=Size, replace=F)	
	tmp <- mclapply(Files[Indx], mc.cores=Cores, mc.preschedule = FALSE, function(X)
	{
		Sample <- ReadSample(X, skip=Skip, SNPList=SNPList, chr=chr)
		if(nrow(Sample) > 1)
		{
			Heterozygosity <- tapply(Sample$B.Allele.Freq, as.factor(Sample$Chr), function(X){ (sum(X > 0.3 & X < 0.7)/length(X))*100 })	
			if(length(Heterozygosity) == length(chr))
			{ 
				return(Heterozygosity)
			} 
		}
	})
	tmp2 <- do.call(rbind.data.frame, tmp)
	colnames(tmp2) <- names(tmp[[1]])
	return(tmp2)
}
