# iPsychCNV

## How to install and test iPsychCNV. 

###In R:
### To install
    setRepositories(ind=1:8)
    # 1 = CRAN , 2 = BioCsoft, 8 = http://R-Forge.R-project.org.
    library(devtools)
    install_github("mbertalan/iPsychCNV")
### Load the package
    library(iPsychCNV)

## Testing iPsychCNV
    # Creating a long mock file (one chromosome).
    LongRoi <- MakeLongMockSample(CNVDistance=1000, Type=c(0,1,2,3,4), Mean=c(-0.6, -0.3, 0.3, 0.6), Size=c(300, 600))
    
    # Running iPsychCNV on long mock data.
    CNVs <- iPsychCNV(PathRawData=".", Pattern="^LongMockSample.tab$", Skip=0)

    # Reading long mock to an object in R.
    Sample <- read.table("LongMockSample.tab", sep="\t", header=TRUE, stringsAsFactors=F)
    
    # Plotting LRR and BAF from 
    PlotLRRAndCNVs(CNVs, Sample, CNVMean=0.3, Name="MyPlotTest.png", Roi=LongRoi)

    # Evaluating CNVs
    CNVs.Eval <- EvaluateMockResults(LongRoi, PredictedCNV)

    # Print ROC curve.
    LongRoc <- plot.roc(CNVs.Eval$CNV.Predicted, CNVs.Eval$CNV.Present, percent=TRUE, print.auc=TRUE)

#### Creating a mock data.
    # Simulates PsychChip array from Illumina.
    #  Creates a Mockfile on local folder and returns the CNV position on the mock sample.
    MockCNVs <­ MockData(N=1, Type="Blood", Cores=1)

#### Predicting CNVs
    CNVs <- iPsychCNV(PathRawData=".", Cores=1, Pattern="*.tab", MINNumSNPs=20, LCR=FALSE, MinLength=10, Skip=0)

#### Subset all CNVs with copy number (CN) different from 2.
    CNVs.Good <- subset(CNVs, CN != 2)
    
#### Creating ROI for Mock Data.
    # ROI: Regions of interest (CNV position on the sample). 
    MockCNVs.Roi <- subset(MockCNVs, CN != 2)
    MockCNVs.Roi$Class <- rep("ROI", nrow(MockCNVs.Roi))

#### Ploting CNVs. 
    # It create a file (test.png) on your local folder.
    PlotAllCNVs(CNVs.Good, Name="test.png", Roi=MockCNVs.Roi)
    
#### Evaluating CNVs
    CNVs.Eval <- EvaluateMockResults(MockCNVs, CNVs.Good)

#### Ploting evaluation using ROC curve.  
    rocobj <- plot.roc(CNVs.Eval$CNV.Predicted, CNVs.Eval$CNV.Present, percent=TRUE,  print.auc=TRUE)  
