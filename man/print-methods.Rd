\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,CNA-method}
\alias{print,DNAcopy-method}
\title{Methods for Function `print' in Package `DNAcopy'}
\description{
  Methods for function \code{print} in Package `DNAcopy'
}
\usage{
\S4method{print}{CNA}(x, \dots)
\S4method{print}{DNAcopy}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"DNAcopy"}}
  \item{...}{For DNAcopy objects there are two optional arguments
    \describe{
      \item{fullOutput:}{logical.  If TRUE the whole output data frame is
	printed and if FALSE only the first 6 rows are.  (default=FALSE)}
      \item{showSegRows:}{logical.  If TRUE the start and end rows of the
	printed segments are shown.  default=FALSE}
    }
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "CNA")}}{
      Prints data.type, number of probes and number of samples.
    }
    \item{\code{signature(x = "DNAcopy")}}{
      Prints the segmentation results.
    }
  }
}

\keyword{methods}
