\name{DNAcopy}
\alias{DNAcopy}
\alias{print.DNAcopy}
\title{Results of segmenting a CNA data object}
\description{
  The results of segmenting data from copy number array experiments from
  programs such as circular binary segmentation (CBS).
}
\value{
  \item{data}{The original CNA object which was the input for segment}
  \item{ID}{sample identifier.}
  \item{chrom}{the chromosome within the sample.}
  \item{loc.start}{the starting map location of the segment}
  \item{loc.end}{the ending map location of the segment}
  \item{num.mark}{the number of markers in the segment}
  \item{data.type}{the segment mean.}
  \item{call}{the call that produced the object.}
}
\details{
  An object of class \code{DNAcopy}.  There is a \code{print} method
  that prints the results in a tabular format.  Each row gives the
  sample, the chromosome, the start and end map locations, the number of
  markers and the mean of each segment.
}
\keyword{nonparametric}
