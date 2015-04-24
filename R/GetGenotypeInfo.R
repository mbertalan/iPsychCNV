GetGenotypeInfo <- function(tmpRaw)
{
	Alleles <- c(tmpRaw$Allele1, tmpRaw$Allele2)
	AllelesPer <- (table(Alleles)/length(Alleles))*100
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
	HeteroPer <- (table(Hetero)/length(Hetero))*100
	
	df <- data.frame(heterozygous=HeteroPer["heterozygous"], homozygous=HeteroPer["homozygous"], GC=GC, AT=AT, "Minus"=AllelesPer["-"], "D"=AllelesPer["D"], "I"=AllelesPer["I"], "A"=AllelesPer["A"], "T"=AllelesPer["T"], "G"=AllelesPer["G"], "C"=AllelesPer["C"], stringsAsFactors=F)
	return(df)
}
