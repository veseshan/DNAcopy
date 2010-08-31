\name{DNAcopy-class}
\docType{class}
\alias{DNAcopy-class}
\alias{$,DNAcopy-method}
\alias{show,DNAcopy-method}
\alias{summary,DNAcopy-method}

\title{Class "DNAcopy"}
\description{The \code{DNAcopy} extends the \code{CNA} class and stores
  the segmentation results in addition to the data.}


\section{Slots}{
  \describe{
    \item{\code{call}:}{the R function call that produced the object.}
    \item{\code{chrom}:}{the chromosome from which the markers come from.
      This variable is numeric or ordered vector.}
    \item{\code{maploc}:}{the locations of the markers along the genome.
      This variable is numeric vector.}
    \item{\code{genomdat}:}{this variable is a matrix, big.matrix or
      ff_matrix.  Each column is a sample and contains the logratio of
      LOH data for each marker. The rows correspond to the markers.
      When genomdat is a big.matrix this variable is an external pointer
      to the data in the backingfile.  For ff_matrix it is an R object
      with the meta data for the matrix in the backingfile.}
    \item{\code{data.type}:}{character variable saying logratio or binary.}
    \item{\code{bigmemorydesc}:}{if genomdat is a matrix or ff_matrix
      this is NULL.  If it is a big.matrix it has the data in the
      description slot of the big.matrix.descriptor.} 
    \item{\code{output}:}{a data frame with the following columns:
      \describe{
	\item{ID:}{sample identifier.}
	\item{chrom:}{the chromosome within the sample.}
	\item{loc.start:}{the starting map location of the segment}
	\item{loc.end:}{the ending map location of the segment}
	\item{num.mark:}{the number of markers in the segment}
	\item{seg.mean:}{the segment mean.}
      }
    }
    \item{\code{segRows}:}{a data frame with start and end rows of segments.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{$}:}{the "$" operator is used to access elements from
      `x' using names.  chrom and maploc, output and segRows can be
      accessed this way.}
    \item{\code{show}:}{\code{signature(x = "DNAcopy")}:
      A one line summary of the object is printed.
    }
  }
}

\author{Venkatraman E. Seshan}
\seealso{\code{\link{segment}}}
\keyword{classes}
