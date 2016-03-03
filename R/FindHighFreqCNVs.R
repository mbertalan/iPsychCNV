##' FindHighFreqCNVs: Identify CNVs with high frequency among the samples. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title FindHighFreqCNVs
##' @param Df: Data frame with predicted CNVs for each sample.
##' @param OverlapCutoff: Unknown, default = 0.8.
##' @return Data frame with predicted CNVs.
##' @author Marcelo Bertalan, Louise K. Hoeffidng. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.

FindHighFreqCNVs <- function(df, OverlapCutoff=0.8)
{
	#tmp <- tapply(df$Start, as.factor(df$Chr), function(X){ sort(table(X), decreasing=TRUE) })
	apply(df, 1, function(X)
	{
		Start1 <- as.numeric(X["Start"])
		Stop1 <- as.numeric(X["Stop"])
		Chr1 <- gsub(" ", "", X["Chr"])
		Length1 <- as.numeric(X["Length"])
		
		apply(df, 1, function(Y)
		{
			Start2 <- as.numeric(Y["Start"])
			Stop2 <- as.numeric(Y["Stop"])
			Chr2 <- gsub(" ", "",Y["Chr"])
			Length2 <- as.numeric(Y["Length"])
	
			if(Chr1 %in% Chr2)
			{
				if(Start1 <= Stop2 & Stop1 >= Start2)	
				{
					if(Start1 >= Start2){ OverStart <- Start1 }else{ OverStart <- Start2 }
					if(Stop1 >= Stop2){ OverStop <- Stop2 }else{ OverStop <- Stop1 }
					
					Overlap <- OverStop - OverStart
					OverL1 <- Overlap/Length1
					OverL2 <- Overlap/Length2
					if(OverL1 > OverlapCutoff & OverL2 > OverlapCutoff)
					{
						return(TRUE)
					}
					else
					{
						return(FALSE)
					}
				}
				else
				{
					return(FALSE)
				}
			}
			else
			{
				return(FALSE)
			}
		})
	})
}	
