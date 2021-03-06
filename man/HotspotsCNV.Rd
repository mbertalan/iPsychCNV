% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HotspotsCNV.R
\name{HotspotsCNV}
\alias{HotspotsCNV}
\title{HotspotsCNV}
\source{
\url{http://biopsych.dk/iPsychCNV}
}
\usage{
HotspotsCNV(df, Freq = 1, OverlapCutoff = 0.7, Cores = 1,
  OverlapMin = 0.9, OverlapMax = 1.1, hg = "hg19")
}
\arguments{
\item{Df:}{Dataframe with CNV predction with chromosome (Chr.), start, and stop position.}

\item{Freq:}{Minimum number of CNVs to be considered a hotspot, default = 1.}

\item{OverlapCutoff:}{Minimum overlap among CNVs to be considered the same CNV region, default = 0.7.}

\item{Cores:}{Numeric, Number of cores used, default = 1.}

\item{OverlapMin:}{Minimum overlap with hotspot to be selected for counting, default = 0.9.}

\item{OverlapMax:}{Maximum overlap with hotspot to be selected for counting, default = 1.1.}

\item{hg:}{Human genome version, default = hg19. For full range of possibilities, run ucscGenomes()$db {rtracklayer}}
}
\value{
CNV hotspots.
}
\description{
HotspotsCNV: Identify Copy Number Variation (CNV) hotspots.
}
\details{
Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
}
\examples{
MockDataCNVs <- MockData(N=100, Type="PKU", Cores=20)
iPsych.Pred <- iPsychCNV(PathRawData=".", MINNumSNPs=20, Cores=20, Pattern="^MockSample", MinLength=10, Skip=0)
iPsych.Pred.hotspots <- HotspotsCNV(iPsych.Pred, Freq=2, OverlapCutoff=0.9, Cores=1)
}
\author{
Marcelo Bertalan, Louise K. Hoeffding, Ida Elken Sønderby.
}

