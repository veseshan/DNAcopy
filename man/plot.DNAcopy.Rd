\name{plot.DNAcopy}
\alias{plot.DNAcopy}
\title{Plot the data and results from segment of a CNA object}
\description{
  Plots the data from a copy number array experiment (aCGH, ROMA etc.)
  along with the results of segmenting it into regions of equal copy
  numbers.
} 
\usage{
  plot.DNAcopy(x, plot.type=c("whole", "samplebychrom", "chrombysample",
            "plateau"),  sbyc.layout = NULL, cbys.nchrom=1, cbys.layout=NULL,
             include.means = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{DNAcopy} resulting from analyzing
    data from copy number array experiments}
  \item{plot.type}{the type of plot.}
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
\keyword{nonparametric}
