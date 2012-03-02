\name{CNA}
\alias{CNA}
\alias{print.CNA}
\title{Create `Copy Number Array' data object}
\description{
  Creates a `copy number array' data object used for DNA copy number
  analyses by programs such as circular binary segmentation (CBS).
}
\usage{
  CNA(genomdat, chrom, maploc, data.type=c("logratio","binary"),
                 sampleid=NULL, presorted = FALSE)
  \method{print}{CNA}(x, \dots)
}
\arguments{
  \item{genomdat}{a vector or matrix of data from array-CGH, ROMA, or
    other copy number experiments. If it is a matrix the rows correspond
    to the markers and the columns to the samples.}
  \item{chrom}{the chromosomes (or other group identifier) from which
    the markers came.  Vector of length same as the number of rows of
    genomdat.  If one wants the chromosomes to be ordered in the
    natural order, this variable should be numeric or ordered category.} 
  \item{maploc}{the locations of marker on the genome.  Vector of length
    same as the number of rows of genomdat. This has to be numeric.}
  \item{data.type}{logratio (aCGH, ROMA, etc.) or binary (LOH).}
  \item{sampleid}{sample identifier.  If missing the samples are named
    by prefixing "Sample" to consecutive integers.}
  \item{presorted}{logical indicator telling if the data have already
    been sorted by chrom and maploc.  Default is FALSE.}
  \item{x}{object returned by CNA}
  \item{...}{arguments to be passed onto print command called within.}
}
\value{
  An object of class \code{CNA}.  There is a \code{print} method that
  gives the number of samples and probes and the type of data.
}

\details{
  Data that are NA, Inf, NaN will be removed on a per sample basis for
  "genomdat" and all samples for "chrom" and "maploc".

  If the chrom variable has non-numeric values make it into an ordered
  variable to get them ordered correctly. E.g. for human genome use:
  \code{chrom <- ordered(chrom, levels=c(1:22,"X","Y"))} to prepare the
  variable if chromosomes X and Y are present in your data.
}

\examples{

data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

}

\keyword{nonparametric}



