##' ClassNumbers: Define how B allele frequency (BAF) is distributed over five diferent classes.
##'
##' Receive a vector with BAF and returns a data frame with a populate percentage for each class. Each class represent a group where a BAF is more likely to be found. 
##' @title ClassNumbers
##' @return Data frame with a populate percentage for each class
##' @author Marcelo Bertalan
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export


ClassNumbers <- function(tmpRaw)
{
	#Centers <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)
	#names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	tmp2 <- subset(tmpRaw, PFB > 0 & PFB < 1)
	UsedBAF <- round((nrow(tmp2)/nrow(tmpRaw))*100, 1)
	Centers <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)
	names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	df <- data.frame(AAAA=0, AAAB=0, AAB=0, AB=0, ABB=0, ABBB=0, BBBB=0)
	if(nrow(tmp2) > 9)
	{
		BAF <- tmp2$B.Allele.Freq
		tmp <- sapply(Centers, function(X){ abs(BAF - X) })
		tmp2 <- apply(tmp, 1, function(X){ names(X)[min(X) == X][1] })
		tmp3 <- data.frame(Class=tmp2, BAF=BAF, Indx=1:length(tmp2), stringsAsFactors=F)
		res <- round((table(tmp3$Class)/sum(table(tmp3$Class)))*100, 1)
		df[names(res)] <- res
	}
	df$UsedBAF <- UsedBAF
	return(df)
}
	
	
	
