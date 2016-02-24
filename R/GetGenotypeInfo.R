##' GetGenotypeInfo: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetGenotypeInfo
##' @param tmpRaw: Unknown.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown
##'

GetGenotypeInfo <- function(tmpRaw)
{
	if(length(tmpRaw$Allele1) > 0 & length(tmpRaw$Allele2) > 0)
	{
		Alleles <- c(tmpRaw$Allele1, tmpRaw$Allele2)
		AllelesPer <- c(0,0,0,0,0,0,0)
		names(AllelesPer) <- c("-","A","C","D","G","I","T")
		AllelesPer2 <- (table(Alleles)/length(Alleles))*100
		AllelesPer[names(AllelesPer2)] <- AllelesPer2
	
		GC <- AllelesPer["G"] +  AllelesPer["C"]
		AT <- AllelesPer["A"] +  AllelesPer["T"]
		Hetero <- apply(tmpRaw[, c("Allele1", "Allele2")], 1, function(X)
		{	
			if(X[1] %in% X[2])
			{
				return("homozygous")
			}	
			else
			{
				return("heterozygous")
			}	
		})
		HeteroPer <- c(0,0)
		names(HeteroPer) <- c("homozygous", "heterozygous")
		HeteroPer2 <- (table(Hetero)/length(Hetero))*100
		HeteroPer[names(HeteroPer2)] <- HeteroPer2
		df <- data.frame(heterozygous=HeteroPer["heterozygous"], homozygous=HeteroPer["homozygous"], GC=GC, AT=AT, "Minus"=AllelesPer["-"], "D"=AllelesPer["D"], "I"=AllelesPer["I"], "A"=AllelesPer["A"], "T"=AllelesPer["T"], "G"=AllelesPer["G"], "C"=AllelesPer["C"], stringsAsFactors=F)
	}
	else
	{
		df <- data.frame(heterozygous=NA, homozygous=NA, GC=NA, AT=NA, "Minus"=NA, "D"=NA, "I"=NA, "A"=NA, "T"=NA, "G"=NA, "C"=NA, stringsAsFactors=F)
	}
	return(df)
}
