##' GetLocus: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title GetLocus
##' @param Df: Data frame with predicted CNVs for each sample.
##' @return data frame with CNV locus.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

GetLocus <- function(df) # df with CHR=chr12, CNV_Start=numeric, CNV_Stop=numeric
{
	library(biovizBase)
	library(BiocGenerics)
	data(hg19IdeogramCyto)
	tmp <- BiocGenerics::as.data.frame(hg19IdeogramCyto)
	df$CHR <- sapply(df$Chr, function(X){ paste("chr", X, sep="", collapse="") })
	if(length(df$CNV_Start) == 0){ 	df$CNV_Start <- df$Start }
	if(length(df$CNV_Stop) == 0){ 	df$CNV_Stop <- df$Stop }
	

	Locus <- apply(df, 1, function(X){ subset(tmp, seqnames %in% X["CHR"] & end > as.numeric(X["CNV_Start"]) & start < as.numeric(X["CNV_Stop"]))$name })
	Locus2 <- lapply(Locus, function(L){ if(length(L) > 1){ paste(droplevels(L), sep="-", collapse="-") }else{ as.character(droplevels(L)) } })
	tmp <- lapply(Locus2, length)
	tmp2 <- unlist(tmp)
	Locus2[tmp2 == 0] <- "NotAvailable"
	Locus2  <- unlist(Locus2)

	df$locus <- Locus2
	Locus <- apply(df, 1, function(X){ paste(X["Chr"],X["locus"], sep="", collapse="") })
	df$locus <- Locus
	return(df)
}
	
