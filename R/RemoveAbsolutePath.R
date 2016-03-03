##' RemoveAbsolutePath: Remove path from Filename with absolute path
##'
##' @title RemoveAbsolutePath
##' @param object - file name with absolute path
##' @return return file name without absolute path
##' @author Ida Elken Sønderby
##' @export
##' @example
##' filename <- "/Volumes/CNVs/BAFLRR_Files/SampleID"
##' RemoveAbsolutePath(filename)

RemoveAbsolutePath <- function(object)
{
  object <- gsub(".+/","",object)
  object <- gsub("/","",object)
  return(object)
}
