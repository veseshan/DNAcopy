\name{subset-methods}
\docType{methods}
\alias{subset-methods}
\alias{subset,CNA-method}
\alias{subset,DNAcopy-method}
\title{Methods for Function `subset' in Package `DNAcopy'}
\description{
  Methods for function \code{subset} in Package `DNAcopy'
}
\usage{
\S4method{subset}{CNA}(x, chromlist=NULL, samplelist=NULL,
     fbdata=NULL, backingfile=NULL, \dots)
\S4method{subset}{DNAcopy}(x, chromlist=NULL, samplelist=NULL,
     fbdata=NULL, backingfile=NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{"CNA"} \code{"DNAcopy"} }
  \item{chromlist}{a vector of chromosome numers or names.  If missing
    the whole genome is returned.}
  \item{samplelist}{the sample for which the plot is requested.  Should
    be a vector of valid sample names or numbers.  If missing all the
    sample are returned.}
  \item{fbdata}{Logical.  If genomdat in x is a big.matrix then this
    options should be specified.  If true then a new file with the
    genomdat for the subset of samples and chromosomes will be created.}
  \item{backingfile}{Charcter.  If genomdat in x is a big.matrix and
    fbdata is TRUE then this option should be specified giving a valid
    filename.  If file exists it will be overwritten.}
  \item{\dots}{other arguments to \code{subset} function. none currently.}    
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "CNA")}}{
      Data subsetted to the chromosomes and samples specified.
    }
    \item{\code{signature(x = "DNAcopy")}}{
      Data subsetted to the chromosomes and samples specified.
    }
  }
}
\keyword{methods}

