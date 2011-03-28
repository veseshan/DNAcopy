\name{exon.segment}
\alias{exon.segment}
\title{Binary segmentation of exon data.}
\description{
  Compute the binary segmentation statistic, location and approximate p-value.
}
\usage{
  exon.segment(gene, eloc, edat, ngrid=100, tol=1e-6)
}
\arguments{
  \item{gene}{gene names in the exon data}
  \item{eloc}{exon locations within gene}
  \item{edat}{exon expressions within gene}
  \item{ngrid}{number grid points for the integral}
  \item{tol}{tolerance level for calculating nu}
}
\value{
  a matrix with three columns.  The maximal statistic from binary
  segmentation, its location and the p-values for each gene.
}
\details{
  The p-values are obtained by applying Siegmund's approximation for the
  maximal statistic from binary segmenting consecutive segments within a
  chromosome.  These are one-sided test for an increase in expression.
}
\examples{

# test code on an easy data set
set.seed(25)
gene <- rep(c("A", "B"), c(30,20))
eloc <- c(1:30, 1:20)
edat <- matrix(rnorm(500), 50, 10)
# changes for gene1 in samples 3 & 7
edat[1:30, 3] <- edat[1:30, 3] + rep(0.9*0:1, c(17, 13))
edat[1:30, 7] <- edat[1:30, 7] + rep(1.1*0:1, c(21, 9))
# changes for gene2 in samples 4 & 7
edat[31:50, 4] <- edat[31:50, 4] + rep(1.1*0:1, c(8, 12))
edat[31:50, 7] <- edat[31:50, 7] + rep(1.2*0:1, c(13, 7))
exon.segment(gene, eloc, edat)

}

\author{Venkatraman E. Seshan}

\keyword{nonparametric}
