\name{CNA-class}
\docType{class}
\alias{CNA-class}
\alias{[<-,CNA,index,index-method}
\alias{[<-,CNA,index,missing-method}
\alias{[<-,CNA,missing,index-method}
\alias{[<-,CNA,missing,missing-method}
\alias{[,CNA,ANY,ANY,ANY-method}
\alias{[,CNA,index,index,missing-method}
\alias{[,CNA,index,missing,missing-method}
\alias{[,CNA,missing,index,missing-method}
\alias{[,CNA,missing,missing,ANY-method}
\alias{$<-,CNA-method}
\alias{$,CNA-method}
\alias{show,CNA-method}

\title{Class "CNA"}
\description{The \code{CNA} class is designed to enable users to handle
  very large filebased data using bigmemory or ff.}
  

\section{Slots}{
\describe{
  \item{\code{chrom}:}{the chromosome from which the markers come from.
    This variable is numeric or ordered vector.}
  \item{\code{maploc}:}{the locations of the markers along the genome.
    This variable is numeric vector.}
  \item{\code{genomdat}:}{this variable is a matrix, big.matrix or
    ff_matrix.  Each column is a sample and contains the logratio of LOH
    data for each marker. The rows correspond to the markers.
    When genomdat is a big.matrix this variable is an external pointer
    to the data in the backingfile.  For ff_matrix it is an R object
    with the meta data for the matrix in the backingfile.}
  \item{\code{data.type}:}{character variable saying logratio or binary.}
  \item{\code{bigmemorydesc}:}{if genomdat is a matrix or ff_matrix this
    is NULL.  If it is a big.matrix it has the data in the description
    slot of the big.matrix.descriptor.}
}}

\section{Methods}{
  \describe{
    \item{\code{$}:}{the "$" operator is used to access elements from
      `x' using names.  chrom and maploc can be accessed this way.}
    \item{\code{$<-}:}{the "$<-" operator is used to replace elements of
      `x' using names.  chrom and maploc can be changed this way.}
    \item{\code{[}:}{the "[" operator is used to extract elements from
      the genomdat matrix.  It takes \code{i} and \code{j} to refer to
      the rows and columns needed.}
    \item{\code{[<-}:}{the "[<-" operator is used to replace elements
      in the genomdat matrix.  It takes \code{i} and \code{j} to refer
      to the rows and columns needed.}
    \item{\code{show}:}{\code{signature(x = "CNA")}:
      A one line summary of the object is printed.
    }
  }
}

\author{Venkatraman E. Seshan}
\seealso{\code{\link{CNA}}}
\keyword{classes}
