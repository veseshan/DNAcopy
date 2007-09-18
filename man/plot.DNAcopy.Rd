\name{plot.DNAcopy}
\alias{plot.DNAcopy}
\title{Plot the data and results from segment of a CNA object}
\description{
  Plots the data from a copy number array experiment (aCGH, ROMA etc.)
  along with the results of segmenting it into regions of equal copy
  numbers.
} 
\usage{
  \method{plot}{DNAcopy}(x, plot.type=c("whole", "plateau", "samplebychrom",
               "chrombysample"), xmaploc=FALSE, altcol=TRUE, sbyc.layout=
                NULL, cbys.nchrom=1, cbys.layout=NULL, include.means=TRUE,
                zeroline=TRUE, pt.pch=NULL, pt.cex=NULL, pt.cols=NULL,
                segcol= NULL, zlcol=NULL, ylim=NULL, lwd=NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{DNAcopy} resulting from analyzing
    data from copy number array experiments.}
  \item{plot.type}{the type of plot.}
  \item{xmaploc}{logical flag to indicate that the X axis is the maploc
    position rather than the index.  Since the segments are rearranged
    the plateau plot does not use maploc position.}
  \item{altcol}{logical flag to indicate if chromosomes should be
    plotted in alternating colors in the whole genome plot.}
  \item{sbyc.layout}{\code{layout} settings for the multifigure grid layout
    for the `samplebychrom' type.  It should be specified as a vector of
    two integers which are the number of rows and columns.  The default
    values are chosen based on the number of chromosomes to produce a
    near square graph.   For normal genome it is 4x6 (24 chromosomes)
    plotted by rows.}
  \item{cbys.layout}{\code{layout} settings for the multifigure grid layout
    for the `chrombysample' type.  As above it should be specified as
    number of rows and columns and the default chosen based on the
    number of samples.}
  \item{cbys.nchrom}{the number of chromosomes per page in the layout.
    The default is 1.}
  \item{include.means}{logical flag to indicate whether segment means
    are to be drawn.}
  \item{zeroline}{logical flag to indicate whether a horizontal line at
    y=0 is to be drawn.}
  \item{pt.pch}{the plotting character used for plotting the log-ratio
    values (default is ".").}
  \item{pt.cex}{the size of plotting character used for the log-ratio
    values (default is 3).}
  \item{pt.cols}{the color list for the points. The colors alternate
    between chromosomes. If missing the point colors are black and green.}
  \item{segcol}{the color of the lines indicating the segment means. If
    missing the line color is set to be red.} 
  \item{zlcol}{the color of the zeroline. If missing it is set to be grey.}
  \item{ylim}{this argument is present to override the default limits
    which is the range of symmetrized log-ratios.}
  \item{lwd}{line weight of lines for segment mean and zeroline.  If
    missing it is set to 3.}
  \item{...}{other arguments which will be passed to \code{plot}
    commands.}
}
\details{
  There are four possible plot types.  For the type `whole' the data
  are plotted for the entire genome.  For the `samplebychrom' type a
  graph with each chromosome (of a given sample) is drawn in a separate
  figure on a multi-figure grid.  For the `plateau' type the graph
  is drawn with the chromosome segments re-ordered by the segment means.
  For the `chrombysample' type the samples for a given chromosome are
  drawn in a 4x6 multi-figure grid in multiples of 24. By default the
  segments means are drawn.  For multisample data each sample or
  chromosome is drawn on a separate sheet.  When invoked interactively
  the user is prompted before advancing to the next sample. 
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

#Plot whole studies

plot(segment.smoothed.CNA.object, plot.type="w")

#Plot each study by chromosome

plot(segment.smoothed.CNA.object, plot.type="s")

#Plot each chromosome across studies (6 per page)

plot(segment.smoothed.CNA.object, plot.type="c", cbys.layout=c(2,1), cbys.nchrom=6)

#Plot by plateaus

plot(segment.smoothed.CNA.object, plot.type="p")

}

\keyword{nonparametric}

