\name{subset.CNA}
\alias{subset.CNA}
\title{Subset a `Copy Number Array' data object}
\description{
   Function to return a subset of a copy number array data object by a
   list of chromosomes and sample.
}
\usage{
  subset.CNA(x, chromlist=NULL, samplelist=NULL, ...)
}                          
\arguments{
  \item{x}{Copy number array data object}
  \item{chromlist}{chromosomes of interest.  Should be a subset of the
    valid chromosome names in the original data.}
  \item{samplelist}{samples of interest.  Can be integers denoting the
    samples of interest or a vector of valid sample names.}
  \item{...}{other arguments which may be passed to \code{subset}.}
}
\value{
  An object of class \code{CNA} with the data for the list of
  chromosomes and samples of interest. 
}
\keyword{nonparametric}
