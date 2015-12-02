##' iPsychCNV: Find Copy Number Variation (CNV) from SNP genotyping arrays. 
##' Specifically designed to reduce false positive CNVs and handle data from amplified DNA on dried blood spots.
##'
##' BAF drift: Calculate the drift in BAF. Used at PeaksAndPits function.
##'
##' @title BAF_drift
##' @param BAF: BAF information from CNV region.
##' @return Data frame with QC variables.
##' @author Marcelo Bertalan
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' BAF <- subCNV$B.Allele.Freq
##' BAF_Res <- BAF_drift(BAF)

BAF_drift <- function(BAF)
{
	# Center
	# AAAA	0.05
	# AAAB	0.3
	# AABB	0.5
	# ABBB	0.7
	# BBBB	0.95
	BAF[is.na(BAF)] <- 0
	Centers <- c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
	names(Centers) <- c("AAAA", "AAAB", "AAB", "AB","ABB", "ABBB", "BBBB") # AB = AABB

	Diff <- sapply(Centers, function(X){  abs(BAF - X) })
	Drift <- apply(Diff, 1, function(X){ c(min(X), names(X)[X == min(X)][1] ) })
	tmp2 <- as.data.frame(t(Drift), stringsAsFactors=F)
	names(tmp2) <- c("Drift", "Class")
	tmp2$Drift <- as.numeric(tmp2$Drift)
	return(tmp2)
}
		
