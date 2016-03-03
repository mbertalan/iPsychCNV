##' RandomlySelectn: Randomly select n from vector or all if length of vector is less than n
##'
##' @title RandomlySelectn
##' @param Individuals - vector
##' @param n - numeric
##' @return n Individuals from vector
##' @author Ida Sønderby
##' @export
##' @example
##' numbers <- seq(1:10)
##' RandomlySelectn(numbers, 2)

# Function to randomly select n individuals but taking into account that if less than ten...
RandomlySelectn <- function(Individuals, n)
{
  if (length(Individuals) <=n)
  {
    return(Individuals)
  }
  else
  {
    return(sample(Individuals,n))
  }
}
