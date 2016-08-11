##' Affymetrix Convert: Calculate LRR and BAF from Affymetix CEL files
##' to allow for CNV detection with iPsych CNV.
##'
##' Compatable with conversion from Affymetrix Assay_6.0 and
##' Assay_5.0. This function applies the instructions given for
##' conversion of Affymetrix data to allow for PennCNV detection from
##' http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/
##'
##' @title Affymetrix Convert
##' @param PathToCel: The path to per sample Affymetrix CEL files.
##' @param Outdir: Output directory for the processed data and conversion files.
##' @param AffyAssay: Character string defining the affymetix assay to be converted (Affy_6.0 or Affy_5.0).
##' @return Data frame with predicted CNVs.
##' @author Johan Hilge Thygesen
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples
##' # Download Affy6.0 example data from hapmap and run analysis
##' system("mkdir -p affy6.0_example")
##' system("wget ftp://anonymous@ftp.ncbi.nlm.nih.gov/hapmap/raw_data/hapmap3_affy6.0/CHEAP.tgz")
##' system("tar xvzf CHEAP.tgz -C affy6.0_example")
##' affyConvert(PathToCel = "affy6.0_example", Outdir = "affyConvert", PathToPennCNV = "/home/mydir/penncnv/", AffyAssay = "Affy_6.0")
##'
##' # Download Affy5.0 example data from ... and run analysis
##' system("mkdir -p affy5.0_example")
##' system("wget http://www.affymetrix.com/support/downloads/data/GW_SNP5_GCOS_CEL_1.zip")
##' system("unzip GW_SNP5_GCOS_CEL_1.zip -d affy5.0_example")
##' affyConvert(PathToCel = "affy5.0_example", Outdir = "affyConvert", PathToPennCNV = "/home/mydir/penncnv/", AffyAssay = "Affy_5.0")

affyConvert <- function(PathToCel, Outdir = "affyConvert", PathToPennCNV, AffyAssay){ 

    ## Variables
    PathToData=paste0(Outdir,"/data") # Data folder where download files will be stored
    AptTool="apt-1.18.2-x86_64-intel-linux"
    
    ## Check system is linux
    if(Sys.info()[1]!="Linux"){
        stop("Non Linux system detected!\nCalculation of LRR and BAF from Affymetrix CEL files must be done on a Linux system")
    }

    ## Generate folders
    system(paste0("mkdir -p ", PathToData)) # if not exists then create folder to store conversion files
    system(paste0("mkdir -p ", Outdir)) # Make output dir
    system(paste0("mkdir -p ", Outdir,"/temp")) # Make temp dir
    system(paste0("mkdir -p ", Outdir,"/converted")) # Make final_out

    ## Generate list of samples to analyse
    samples <- data.frame(cel_files = list.files(PathToCel, pattern = ".CEL", full.names = T))
    write.table(samples, paste0(PathToCel,"/sample_list.txt"), row.names = F, quote = F)

    ## Download data files needed for conversion if not allready present
    ## PennCNV gw6 pack
    if(!file.exists(paste0(PathToData,"/gw6"))){
        print(" ---- DOWNLADING gw6.tar.gz ---- ")
        system(paste0("wget http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz -P ", PathToData))
        system(paste0("tar -xvzf ",PathToData, "/gw6.tar.gz -C ", PathToData))
    }

    ## Download APT tools and assay specific library files
    if(!AffyAssay %in% c("Affy_6.0", "Affy_5.0")){
        stop("Please indicate a supported Affymetrix assay Affy_6.0 or Affy_5.0")
    }

    ## Affy_6.0
    if(AffyAssay=="Affy_6.0"){
        affy6.0AptFiles <- c(paste0(PathToData,"/",AptTool,".zip"),paste0(PathToData, "/genomewidesnp6_libraryfile.zip"))
        if(!all(file.exists(affy6.0AptFiles))){
            stop(paste0("\n------------------------------------------------------------------------------
  Please download the Affymetrix Power Tools (APT) software package and library files from:

    http://www.affymetrix.com/partners_programs/programs/developer/apt_download/apt_terms.affx?v=apt1.18.2_linux64
    http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip

  YOU WILL NEED to log into the website to download the software (registration is free).

  Copy the downloaded files (", AptTool, " & genomewidesnp6_libraryfile.zip) to ", PathToData,"\n\n"))
        }else if(!all(file.exists(c(paste0(PathToData, "/", AptTool), paste0(PathToData, "/CD_GenomeWideSNP_6_rev3"))))){ # Unpack if not allready done
            print(" ---- UNZIPPING Affymetrix Power Tools & Affy 6.0 library files ---- ")
            system(paste0("unzip ",PathToData, "/", AptTool, " -d ", PathToData))
            system(paste0("unzip ",PathToData, "/genomewidesnp6_libraryfile.zip -d ", PathToData))
        }

        ## Run conversion steps
        ## Substep 1.1 - Generate genotyping calls from CEL files (Min mem is 8)
        system(paste0(PathToData, "/", AptTool, "/bin/apt-probeset-genotype -c ", PathToData, "/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf",
               " -a birdseed --read-models-birdseed ", PathToData, "/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.birdseed.models",
               " --special-snps ", PathToData, "/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.specialSNPs",
               " --out-dir ", Outdir, "/temp",
               " --cel-files ", PathToCel, "/sample_list.txt"))

        ## Substep 1.2 - Allele-specific signal extraction from CEL files
        system(paste0(PathToData, "/", AptTool, "/bin/apt-probeset-summarize --cdf-file ",  PathToData, "/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf",
               " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true",
               " --target-sketch ", PathToData, "/gw6/lib/hapmap.quant-norm.normalization-target.txt",
               " --out-dir ", Outdir, "/temp",
               " --cel-files ", PathToCel, "/sample_list.txt"))

        ## Substep 1.3 - Generate canonical genotype clustering file
        system(paste0("fgrep male ",  Outdir, "/temp/birdseed.report.txt | cut -f 1,2 > ", Outdir, "/gender_file.txt"))
        system(paste0(PathToData, "/gw6/bin/generate_affy_geno_cluster.pl ", Outdir, "/temp/birdseed.calls.txt ", Outdir, "/temp/birdseed.confidences.txt ", Outdir, "/temp/quant-norm.pm-only.med-polish.expr.summary.txt",
               " --locfile ", PathToData, "/gw6/lib/affygw6.hg18.pfb",
               " --sexfile ", Outdir, "/gender_file.txt",
               " --out ", Outdir, "/temp/gw6.genocluster"))

        ## Substep 1.4 - LRR and BAF calculations
        system(paste0(PathToData, "/gw6/bin/normalize_affy_geno_cluster.pl ", Outdir, "/temp/gw6.genocluster ", Outdir, "/temp/quant-norm.pm-only.med-polish.expr.summary.txt",
               " --locfile ", PathToData, "/gw6/lib/affygw6.hg18.pfb",
               " --out ", Outdir, "/temp/gw6.lrr_baf.txt"))

        ## Substep 2.0 Split gw6.lrr_baf.txt with penncnv/kcolumn.pl
        system(paste0(PathToPennCNV, "/kcolumn.pl ", Outdir, "/temp/gw6.lrr_baf.txt split 2 -tab -head 3 -name --beforestring .CEL -start 1 -end 500 -out ", Outdir, "/converted/out"))
    }

    ## Affy_5.0
    if(AffyAssay=="Affy_5.0"){
        affy5.0AptFiles <- c(paste0(PathToData,"/",AptTool,".zip"), paste0(PathToData, "/genomewidesnp5_libraryfile_rev1.zip"), paste0(PathToData, "/GenomeWideSNP_5.r2.zip"))
        if(!all(file.exists(affy5.0AptFiles))){
            stop(paste0("\n------------------------------------------------------------------------------
  Please download the Affymetrix Power Tools (APT) software package and library files from:

    http://www.affymetrix.com/partners_programs/programs/developer/apt_download/apt_terms.affx?v=apt1.18.2_linux64
    http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp5_libraryfile_rev1.zip
    http://www.affymetrix.com/Auth/support/downloads/library_files/GenomeWideSNP_5.r2.zip

  YOU WILL NEED to log into the website to download the software (registration is free).

  Copy the downloaded files (", AptTool, " & genomewidesnp6_libraryfile.zip) to ", PathToData,"\n\n"))
        }else if(!all(file.exists(c(paste0(PathToData, "/", AptTool), paste0(PathToData, "/CD_GenomeWideSNP_5_rev1"), paste0(PathToData, "/GenomeWideSNP_5.r2"))))){ # Unpack if not allready done
            print(" ---- UNZIPPING Affymetrix Power Tools & Affy 5.0 library files ---- ")
            system(paste0("unzip ",PathToData, "/", AptTool, " -d ", PathToData))
            system(paste0("unzip ",PathToData, "/genomewidesnp5_libraryfile_rev1.zip -d ", PathToData))
            system(paste0("unzip ",PathToData, "/GenomeWideSNP_5.r2.zip -d ", PathToData))
        }

        ## Run conversion steps
        ## Substep 1.1 - Generate genotyping calls from CEL files (Min mem is 8)
        system(paste0(PathToData, "/", AptTool, "/bin/apt-probeset-genotype -c ", PathToData, "/GenomeWideSNP_5.r2/GenomeWideSNP_5.Full.r2.cdf",
                      " --chrX-snps ", PathToData, "/CD_GenomeWideSNP_5_rev1/Full/GenomeWideSNP_5/LibFiles/GenomeWideSNP_5.Full.chrx",
                      " --read-models-brlmmp ", PathToData, "/CD_GenomeWideSNP_5_rev1/Full/GenomeWideSNP_5/LibFiles/GenomeWideSNP_5.models",
                      " -a brlmm-p",
                      " --out-dir ", Outdir, "/temp",
                      " --cel-files ", PathToCel, "/sample_list.txt"))
                
        ## Substep 1.2 - Allele-specific signal extraction from CEL files
        system(paste0(PathToData, "/", AptTool, "/bin/apt-probeset-summarize --cdf-file ", PathToData, "/GenomeWideSNP_5.r2/GenomeWideSNP_5.Full.r2.cdf",
               " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true",
               " --target-sketch ", PathToData, "/gw6/libgw5/agre.quant-norm.normalization-target.txt",
               " --out-dir ", Outdir, "/temp",
               " --cel-files ", PathToCel, "/sample_list.txt"))
        
        ## Substep 1.3 - Generate canonical genotype clustering file
        system(paste0("fgrep male ",  Outdir, "/temp/birdseed.report.txt | cut -f 1,2 > ", Outdir, "/gender_file.txt"))
        system(paste0(PathToData, "/gw6/bin/generate_affy_geno_cluster.pl ", Outdir, "/temp/brlmm-p.calls.txt ", Outdir, "/temp/brlmm-p.confidences.txt ", Outdir, "/temp/quant-norm.pm-only.med-polish.expr.summary.txt",
               " --locfile ", PathToData, "/gw6/libgw5/affygw5.hg18.pfb",
               " --sexfile ", Outdir, "/gender_file.txt",
               " --out ", Outdir, "/temp/gw5.genocluster"))

        ## Substep 1.4 - LRR and BAF calculations
        system(paste0(PathToData, "/gw6/bin/normalize_affy_geno_cluster.pl ", Outdir, "/temp/gw5.genocluster ", Outdir, "/temp/quant-norm.pm-only.med-polish.expr.summary.txt",
               " --locfile ", PathToData, "/gw6/libgw5/affygw5.hg18.pfb",
               " --out ", Outdir, "/temp/gw5.lrr_baf.txt"))

        ## Substep 2.0 Split gw6.lrr_baf.txt with penncnv/kcolumn.pl
        system(paste0(PathToPennCNV, "/kcolumn.pl ", Outdir, "/temp/gw5.lrr_baf.txt split 2 -tab -head 3 -name --beforestring .CEL -start 1 -end 500 -out ", Outdir, "/converted/out"))

    }

    ## WARNINGS
    if(nrow(samples)<100){
        cat("\nWARNING: You have fewer than 100 samples!\n For this protocol to work, one need to use at least 100 CEL files to generate a reasonably good clustering file\n see http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/\n for further instructions on small sample sizes\n")
    }
    cat("\nWARNING: This analysis have been done with standard settings which works well in moste cases. \n OBS! your sample might need special treatment. \n Inspect the output of the commands and fix pending issues before proceding!!\n\n")

}


