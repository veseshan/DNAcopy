\name{plotSample}
\alias{plotSample}
\title{Plot the data and results from segmentation for a single sample}
\description{
  Plots the data for a single sample from a copy number array experiment
  (aCGH, ROMA etc.) along with the results of segmenting it into regions
  of equal copy numbers.
} 
\usage{
  plotSample(x, sampleid=NULL, chromlist=NULL, xmaploc=FALSE,
             col=c("black","green"), pch=".", cex=NULL, altcol=TRUE,
             segcol="red", lwd=3, zeroline=TRUE, zlcol="grey",
             xlab=NULL, ylab=NULL, main=NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{DNAcopy} resulting from analyzing
    data from copy number array experiments.}
  \item{sampleid}{the sample for which the plot is requested.  Should be
    a valid sample name or number.  If missing the first sample is plotted.} 
  \item{chromlist}{a vector of chromosome numers or names to be plotted.
    If missing the whole genome is plotted.}
  \item{xmaploc}{a logical indicating if data are plotted against genomic
    position or Index.  Defaults to \code{FALSE}.}
  \item{col}{a vector of two colors that can be used for alternating
    colors for successive chromosomes.}
  \item{pch}{the plotting character.  Defaults to \code{.}.}
  \item{cex}{the size of plotting character.  If missing it is set to 3
    if pch is `.' and 1 otherwise.}
  \item{altcol}{a logical indicating if colors of successive chromosomes
    should be alternated.  Defaults to \code{TRUE}.}
  \item{segcol}{color for segment means.}
  \item{zeroline}{a logical indicating if the zeroline is drawn.
    Defaults to \code{TRUE}.}
  \item{zlcol}{color for zero line.}
  \item{lwd}{thickness of the lines.}
  \item{xlab}{the x-axis lavel. If missing Index or Genomic Position
    will be used depending on xmaploc.}
  \item{ylab}{the y-axis label. If missing log(CN) or LOH will be used
    depending on data type.}
  \item{main}{the main title.  If missing sample name will be used.}
  \item{...}{other arguments to the \code{plot} function can be passed here.}
}
\details{
  This function plots the whole genome and segmentation results for a single
  sample.  This function overcomes the deficiency in the plot.DNAcopy function
  which cycles through all the samples.  If sampleid is not specified the
  first sample is plotted.
}

\examples{

#Read in two examples from Snijders et al.

data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

#We generally recommend smoothing single point outliers before analysis
#Make sure to check that the smoothing is proper

smoothed.CNA.object <- smooth.CNA(CNA.object)

#Segmentation at default parameters

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

# Plot whole sample c13330

plotSample(segment.smoothed.CNA.object, sampleid="c13330")

# Plot only chromosomes 1,3,5,7,9 from first sample
plotSample(segment.smoothed.CNA.object, sampleid=1, chromlist=c(1,3,5,7,9))

}

\keyword{nonparametric}

