% modeled after stats4/man/plot-methods.Rd
\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,DNAcopy-method}
\title{Methods for Function `plot' in Package `DNAcopy'}
\description{
  Plot segmentation results for \code{"DNAcopy"} objects.
}
\usage{
\S4method{plot}{DNAcopy}(x, sampleid=NULL, chromlist=NULL, xmaploc=FALSE,
     col=c("black","green"), pch=".", cex=1, altcol=TRUE, segcol="red",
     lwd=3, zeroline=TRUE, zlcol="grey", xlab="", ylab="", main="", \dots)
}
\arguments{
  \item{x}{an object of class \code{"DNAcopy"} }
  \item{sampleid}{the sample for which the plot is requested.  Should be
    a valid sample name or number.  If missing the first sample is plotted.} 
  \item{chromlist}{a vector of chromosome numers or names to be plotted.
    If missing the whole genome is plotted.}
  \item{xmaploc}{a logical indicating if data are plotted against genomic
    position or Index.  Defaults to \code{FALSE}.}
  \item{col}{a vector of two colors that can be used for alternating
    colors for successive chromosomes.}
  \item{pch}{the plotting character.  Defaults to \code{.}.}
  \item{cex}{the size of plotting character.  Default is 1 (3 of `.').}
  \item{altcol}{a logical indicating if colors of successive chromosomes
    should be alternated.  Defaults to \code{TRUE}.}
  \item{segcol}{color for segment means.}
  \item{zeroline}{a logical indicating if the zeroline is drawn.
    Defaults to \code{TRUE}.}
  \item{zlcol}{color for zero line.}
  \item{lwd}{thickness of the lines.}
  \item{xlab}{the x-axis lavel.}
  \item{ylab}{the y-axis label.}
  \item{main}{the main title.  Default is the sample name.}
  \item{\dots}{other arguments to the \code{plot} function can be passed here.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "DNAcopy")}}{Plot segmented
      results for a single sample and a subset of chromosomes in \code{x}.}
  }
}
\keyword{methods}
