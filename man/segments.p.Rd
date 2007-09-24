\name{segments.p}
\alias{segments.p}
\title{p-values for the change-points}
\description{
  This program computes a p-value for the change-points found by the 
  circular binary segmentation (CBS) algorithm.
}
\usage{
  segments.p(x, ngrid=100, tol=1e-6)
}
\arguments{
  \item{x}{an object of class DNAcopy}
  \item{ngrid}{number grid points for the integral}
  \item{tol}{tolerance level for calculating nu}
}
\value{
  a data frame with eight columns.  The seventh column is the maximal
  statistic from binary segmentation and eighth is the p-values added to
  the six columns from the segment command.
}
\details{
  The p-values are obtained by applying Siegmund's approximation for the
  maximal statistic from binary segmenting consecutive segments within a
  chromosome.  This p-value is only to give the relative importance of
  the change-points as the CBS is different from the algorithm used here.
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

\author{Venkatraman E. Seshan \email{olshena@mskcc.org} }

\keyword{nonparametric}
