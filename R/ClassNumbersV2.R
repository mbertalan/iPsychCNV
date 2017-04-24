ClassNumbersV2 <- function(tmpRaw)
{
	Centers <- c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1)
	names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	df <- data.frame(AAAA=0, AAAB=0, AAB=0, AB=0, ABB=0, ABBB=0, BBBB=0)

	BAF <- tmpRaw$B.Allele.Freq
	Size <- length(BAF) 
	df$AAAA <- round((sum(BAF < 0.1)/Size)*100)
	df$AAAB <- round((sum(BAF > 0.1 & BAF < 0.2)/Size)*100)
	df$AAB <- round((sum(BAF > 0.2 & BAF < 0.45)/Size)*100)
	df$AB <- round((sum(BAF > 0.45 & BAF < 0.55)/Size)*100)
	df$ABB <- round((sum(BAF > 0.55 & BAF < 0.8)/Size)*100)
	df$ABBB <- round((sum(BAF > 0.8 & BAF < 0.9)/Size)*100)
	df$BBBB <- round((sum(BAF > 0.9)/Size)*100)

	return(df)
}
