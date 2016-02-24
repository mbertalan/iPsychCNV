##' CleaningPeaks: Remove noise from peaks in B Allele Frequency (BAF). 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title CleaningPeaks
##' @param Tp: Unknown.
##' @return Peaks for BAF.
##' @author Marcelo Bertalan, Louise K. Hoeffding. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.

CleaningPeaks <- function(tp)
{
	tp$peaks[tp$points < 0.15] <- FALSE
	A <- sum(tp$peaks[tp$pos < 101])
	if(A > 1)
	{
		tmpIndxPeak <- sort(tp$points[tp$pos < 100][tp$peaks[tp$pos < 100]], decreasing=T)[1] # Get the highest peak for range 0-100.
		Pos <- tp$pos[tp$pos < 100][tp$points[tp$pos < 100] == tmpIndxPeak]
		tp$peaks[tp$pos < 100] <- FALSE
		tp$peaks[tp$pos == Pos] <- TRUE
	}
	A <- sum(tp$peaks[tp$pos > 399])
	if(A > 1)
	{
		tmpIndxPeak <- sort(tp$points[tp$pos > 399][tp$peaks[tp$pos > 399]], decreasing=T)[1] # Get the highest peak for range 399 - 500.
		Pos <- tp$pos[tp$pos > 399][tp$points[tp$pos > 399] == tmpIndxPeak]
		tp$peaks[tp$pos > 399] <- FALSE
		tp$peaks[tp$pos == Pos] <- TRUE
	}
	
	A <- sum(tp$peaks)
	if(A == 3)
	{
		B <- sum(tp$peaks[tp$pos > 210 & tp$pos < 290]) # Check if the peak is in the middle.
		if(B == 0) # if not, clean peak that is noise.
		{
			tp$peaks[tp$pos > 100 & tp$pos < 210] <- FALSE
			tp$peaks[tp$pos > 290 & tp$pos < 400] <- FALSE
		}
	}

	return(tp)
}
