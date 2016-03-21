##' ReadSample: Read data. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title ReadSample
##' @param RawFile: The raw data files, default = Test.txt.
##' @param Skip: Integer, the number of lines of the data file to skip before beginning to read data, default = 0.
##' @param LCR: List low copy repeat regions, list of SNPs that should be removed, default = FALSE. 
##' @param PFB: Vector, population frequency 0 to 1 for each SNP in the array, default = NULL.
##' @param Chr: Character, select a specific chromosome to be analyzed, default = NA.
##' @param SNPList: Getting Chr. and Position from another source than the RawFile - input should be the full path of the SNPList with columns: Name, Chr, amd Position. Any positions from the RawFile will be erased. A PFB-column is also allowed but will be overwritten by the PFB-parameter or exchanged with 0.5  
##' @return Return data in a data frame.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown.
##'

ReadSample <- function(RawFile="Test.txt", skip=0, LCR=FALSE, PFB=NULL, chr=NA, SNPList=NULL)
{
  suppressPackageStartupMessages(library(data.table))
  Sample <- fread(RawFile, head=T, sep="\t", skip=skip, verbose=FALSE)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- gsub(" ", ".", colnames(Sample))
  colnames(Sample)[colnames(Sample) %in% "Name"] <- "SNP.Name"
  colnames(Sample)[colnames(Sample) %in% "Chromosome"] <- "Chr"
  colnames(Sample)[colnames(Sample) %in% "Allele1.-.Top"] <- "Allele1"
  colnames(Sample)[colnames(Sample) %in% "Allele2.-.Top"] <- "Allele2"
  colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
  colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
  colnames(Sample)[grep("Log.R.Ratio", colnames(Sample))] <- "Log.R.Ratio" # Remove text from PennCNV format.
  colnames(Sample)[grep("B.Allele.Freq", colnames(Sample))] <- "B.Allele.Freq"
  colnames(Sample)[grep("Chr", colnames(Sample))] <- "Chr"
  colnames(Sample)[grep("Position", colnames(Sample))] <- "Position"
  #CNV <- CNV[,c("SNP.Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "Allele1", "Allele2")] # SNP.Name
  
  # Genotype together (deCODE)
  #if(!is.null(Sample$Genotype))
  #{
  #	Allele <- sapply(Sample$Genotype, function(X){ data.frame(Allele1= unlist(strsplit(X, ""))[1], Allele2=unlist(strsplit(X, ""))[2], stringsAsFactors=F) })
  #	Allele <- MatrixOrList2df(Allele)
  #	Sample$Allele1 <- Allele$Allele1
  #	Sample$Allele2 <- Allele$Allele2
  #}

  # SNP-positionFile
  if(!is.null(SNPList)){
    suppressPackageStartupMessages(library(data.table))
    SNPlist <- fread(SNPList, head=T, skip=skip, verbose=FALSE, stringsAsFactors=FALSE, colClasses =c("Chr"="factor"))
    SNPlist <- as.data.frame(SNPlist)
    colnames(SNPlist)[colnames(SNPlist) %in% "Name"] <- "SNP.Name"
    colnames(SNPlist)[colnames(SNPlist) %in% "Chromosome"] <- "Chr" 
    colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
    colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
    Sample$Chr <- NULL # remove previous values of Chr if present
    Sample$Position <- NULL # remove previous values of Position if present
    total <- merge(Sample, SNPlist,by="SNP.Name")
    Sample <- total    
   }

  # removing chr from chromosome name (deCODE)
  Sample$Chr <- gsub("chr", "", Sample$Chr)
  
  
  # PFB
  if(is.null(PFB)){ Sample$PFB <- rep(0.5, nrow(Sample)) }else{ Sample$PFB <- PFB }
  
  # Subsetting 
  Sample <- subset(Sample, !Chr %in% c("XY", "0")) # "MT", "X", "Y",
  Sample <- subset(Sample, !is.na(Sample$B.Allele.Freq)) # Removes SNPs without BAF-value 
  Sample <- subset(Sample, !is.na(Sample$Log.R.Ratio)) #  Removes SNPs without LRR-value 
  
  # chr specific. Example chr="22"
  if(!is.na(chr)){ Sample <- subset(Sample, Chr %in% chr) }
  
  if(LCR == TRUE)
  {
    Sample <- subset(Sample, !SNP.Name %in% LCR.SNPs) 
  }
  
  Sample$LRR <- Sample$Log.R.Ratio # CNV$LRR is the
  return(Sample)
}
