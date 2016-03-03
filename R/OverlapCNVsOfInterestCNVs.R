##' OverlapCNVOfInterestCNVs: Find overlap between a list of CNVs of interest and a query with CNVs
##'
##' @title OverlapCNVOfInterestCNVs
##' @param cnvsofinterest - Dataframe of CNVsofInterest with as minimum columns Chr, start, stop [Chr should be coded with numbers]
##' @param query - Dataframe of CNVs with as minimum columns Chr, Start, Stop, Length & CN [Chr should be coded with numbers]
##' @return return merged list of cnvsofinterest and query if there is a match/overlap.
##' @author Ida Sønderby
##' @export
##' @example
##' cnvsofinterest <- data.frame(Locus = c("Chr1_110000-120000", "Chr1_120000-170000"), Chr = c(1,1), start = c(110000, 120000), stop=c(120000, 170000))
##' cnvs <- data.frame(Chr = c(1,1), Start = c(100000, 130000), Stop = c(120000, 160000), Length = c(20000, 30000), CN = c(1,1) )
##' OverlapCNVOfInterestCNVs(cnvsofinterest, cnvs)

OverlapCNVOfInterestCNVs <- function(cnvsofinterest, query)
{
  # obtain package required
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  require(GenomicRanges)

  # make a query group
  gr1 = with(cnvsofinterest, GRanges(Rle(factor(Chr,
                                                levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))),
                                     IRanges(start, stop)))

  # make a subjectgroup
  gr2 = with(query, GRanges(Rle(factor(Chr,
                                       levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))),
                            ranges = IRanges(Start, Stop)))

  # find overlap between query and subject
  olaps = findOverlaps(gr1, gr2)

  # couple a list of query and subject
  result = data.frame(cnvsofinterest[queryHits(olaps),], query[subjectHits(olaps),])
  return(result)
}
