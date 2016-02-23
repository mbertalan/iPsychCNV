##' ConvertPenntoRCNVs: converts PennCNV-output to R-readable file with the extension ".Rtxt" and reads in as dataframe.
##'
##' Specifically designed to handle noisy data from amplified DNA on Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ConvertPenntoRCNVs
##' @param PennCNVFile: filepath with CNVs called in PennCNV
##' @return Return output from PennCNV in dataframe and a R-readable file with the extension ".Rtxt".
##' @author Ida Sønderby, Marcelo Bertalan, Louise K. Hoeffding
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples unknown
##' 


ConvertPenntoRCNVs <- function(PennCNVFile = "Test_CNVs.txt") {

  # Load conversion script (pennCNV-format to R-format)
  Penn2Tab <- system.file("exec/Penn2Tab.pl",package="iPsychCNV")


  # Make a new name for the output-data
  CNVs_R <- paste(PennCNVFile, ".Rtxt", sep="", collapse="") # extension=".Rtxt"

  # convert data-command
  Command <- paste(Penn2Tab, " < ", PennCNVFile, " > ", CNVs_R, sep="", collapse="")
  cat(Command, "\n")

  # run command that creates a file in the same folder as filtered QCed CNVs with addition ".Rtxt"
  system(Command)

  CNVs <- read.table(CNVs_R, sep="\t", header=TRUE) # read in data in dataframe
  CNVs$Source <- as.factor(rep("PennCNV", nrow(CNVs))) # add source of CNVs as PennCNVs
  CNVs$ID <- CNVs$File # adds an ID-column
  return(CNVs)
}
