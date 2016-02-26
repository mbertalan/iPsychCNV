##' format_si: Format a vector of numeric values according to the International System of Units.
##'
##' Format a vector of numeric values according to the International System of Units.
##' http://en.wikipedia.org/wiki/SI_prefix
##' Based on code by Ben Tupper
##' https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
##' Args:
##'   ...: Args passed to format()
##'
##' @title format_si
##' @param x: A vector of numeric values
##' @return A function to format a vector of strings using SI prefix notation
##' @author Stefan Möding (http://www.moeding.net/archives/32-Metric-prefixes-for-ggplot2-scales.html)
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown

format_si <- function(...) {
  function(x) {
    limits <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12,
                1e-9,  1e-6,  1e-3,  1e0,   1e3,
                1e6,   1e9,   1e12,  1e15,  1e18,
                1e21,  1e24)
    prefix <- c("y",   "z",   "a",   "f",   "p",
                "n",   "µ",   "m",   " ",   "k",
                "M",   "G",   "T",   "P",   "E",
                "Z",   "Y")

    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)

    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)

    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...),
          prefix[i])
  }
}
