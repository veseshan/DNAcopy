\name{segments.p}
\alias{segments.p}
\title{p-values for the change-points}
\description{
  This program computes pseudo p-values and confidence intervals for the
  change-points found by the circular binary segmentation (CBS) algorithm.
}
\usage{
  segments.p(x, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)
}
\arguments{
  \item{x}{an object of class DNAcopy}
  \item{ngrid}{number grid points for the integral}
  \item{tol}{tolerance level for calculating nu}
  \item{alpha}{Confidence level is 1-alpha}
  \item{search.range}{statistic is maximized over nu +/- search.range}
  \item{nperm}{number of permutations for confidence interval} 
}
\value{
  a data frame with ten columns.  The maximal statistic from binary
  segmentation, the p-values and lower and upper alpha/2 confidence
  limits (as genomic positions) are added to the six columns from the
  segment command.

  NOTE:  THE p VALUES ARE APPROXIMATE TAIL PROBABILITIES.  ANY VALUE
  GREATER THAN 0.1 CAN HAVE LARGE ERROR.  p > 1 ARE REPLACED WITH 1.
}
\details{
  The p-values are obtained by applying Siegmund's approximation for the
  maximal statistic from binary segmenting consecutive segments within a
  chromosome.  This p-value is only to give the relative importance of
  the change-points as the CBS is different from the algorithm used here.

  The confidence intervals are obtained by a permutation algorithm.  The
  data are permuted to the left and right of the identified change-point
  and the location of the maximal binary segmentation statistic computed.
  The confidence interval is given by the quantiles of the permutation
  distribution of the locations.

  The statistical properties of this confidence interval is unknown.  It
  is used to give an idea of the uncertainity on the location of the
  change-points as the CBS is different from the algorithm used here.
}
\examples{

# test code on an easy data set
set.seed(25)
genomdat <- rnorm(500, sd=0.1) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
plot(genomdat)
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test1 <- segment(CNA(genomdat, chrom, maploc))
segments.p(test1)

}

\author{Venkatraman E. Seshan}

\keyword{nonparametric}
