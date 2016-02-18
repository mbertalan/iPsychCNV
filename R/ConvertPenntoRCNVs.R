##' ConvertPenntoRCNVs: converts PennCNV-output to R-readable file with the extension ".Rtxt" and reads in as dataframe
##'
##' Specifically designed to handle noisy data from amplified DNA on  Phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ConvertPenntoRCNVs
##' @return return output from PennCNV in dataframe
##' @param PennCNVFile - file with CNVs called in PennCNV
##' @author Marcelo Bertalan/Ida Sønderby
##' @export

# Function that takes pennCNV-output and convert to R-readable file in the same folder with extension ".Rtxt" and reads in this dataframe
ConvertPenntoRCNVs <- function(File = PennCNVFile) { 
  
  # load Conversion script (pennCNV-format to R-format)
  Penn2Tab <- system.file("exec/Penn2Tab.pl",package="iPsychCNV")
  
  
  # Make a new name for the output-data
  CNVs_R <- paste(PennCNVFile, ".Rtxt", sep="", collapse="") # a way to get output in particular place...
  
  # convert data-command
  Command <- paste(Penn2Tab, " < ", PennCNVFile, " > ", CNVs_R, sep="", collapse="")
  cat(Command, "\n")
  
  # run command that creates a file in the same folder as filtered QCed CNVs with addition .Rtxt
  system(Command) 
  
  CNVs <- read.table(CNVs_R, sep="\t", header=TRUE) # read in data
  CNVs$Source <- as.factor(rep("PennCNV", nrow(CNVs))) # add source of CNVs
  CNVs$ID <- CNVs$File
  return(CNVs)
}
