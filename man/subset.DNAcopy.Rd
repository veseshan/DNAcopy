\name{subset.DNAcopy}
\alias{subset.DNAcopy}
\title{Subset a DNAcopy data object}
\description{
   Function to return a subset of a copy number array data object by a
   list of chromosomes and sample.
}
\usage{
  \method{subset}{DNAcopy}(x, chromlist=NULL, samplelist=NULL, ...)
}                          
\arguments{
  \item{x}{DNAcopy object}
  \item{chromlist}{chromosomes of interest.  Should be a subset of the
    valid chromosome names in the original data.}
  \item{samplelist}{samples of interest.  Can be integers denoting the
    samples of interest or a vector of valid sample names.}
  \item{...}{other arguments which may be passed to \code{subset}.}
}
\value{
  An object of class \code{DNAcopy} with the input data and the results
  of segmenting them only for the chromosomes and samples of interest.
}
\keyword{nonparametric}
