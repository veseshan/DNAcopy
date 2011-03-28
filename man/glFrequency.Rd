\name{glFrequency}
\alias{glFrequency}
\title{Additional summary measured for the segments}
\description{
  This program computes the frequency of gains and losses for each probe
  as a function of level of mad.
}
\usage{
  glFrequency(xout, threshold=1)
}
\arguments{
  \item{xout}{an object of class DNAcopy}
  \item{threshold}{threshold value to call gain or loss}
}
\value{
  A segment is called a gain or loss if the segment mean is at least the
  threshold* mad distance away from the median copy number level.  The
  output is a data frame with five columns which give the chromosome
  (chrom), genomic position (maploc), the number of samples with
  available data (pfreq), and the gain (gain) and loss (loss).
}

\author{Venkatraman E. Seshan}

\keyword{nonparametric}
