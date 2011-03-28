\name{segments.summary}
\alias{segments.summary}
\title{Additional summary measured for the segments}
\description{
  This program computes the standard deviation, median and the mad of
  the data for each segment found by the CBS algorithm.
}
\usage{
  segments.summary(x)
}
\arguments{
  \item{x}{an object of class DNAcopy}
}
\value{
  a data frame with nine columns.  The sd, median and mad of each
  segment is added to the six columns from the segment command.
}
\examples{

# test code on an easy data set
set.seed(25)
genomdat1 <- rnorm(500, sd=0.1) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
genomdat2 <- rnorm(500, sd=0.1) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
genomdat1[sample(1:500,5)] <- NA
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test1 <- segment(CNA(cbind(genomdat1,genomdat2), chrom, maploc))
segments.summary(test1)

}

\author{Venkatraman E. Seshan}

\keyword{nonparametric}
