\name{segment}
\alias{segment}
\title{Genome Segmentation Program}
\description{
  This program segments DNA copy number data into regions of estimated 
  equal copy number using circular binary segmentation (CBS).
}
\usage{
  segment(genomdat, chrom, maploc, data.type = c("logratio",
                 "binary"), alpha = 0.01, nperm = 10000, window.size =
                 NULL, overlap = 0.25, trim = 0.025, smooth.outliers =
                 TRUE, smooth.region = 2, outlier.SD = 4, smooth.SD =
                 2, smooth.output = FALSE, undo.splits = c("none",
                 "prune", "sdundo"), undo.prune = 0.05, undo.SD = 3,
                 verbose = TRUE)
}
\arguments{
  \item{genomdat}{a vector or matrix of data from array-CGH, ROMA, or
    other copy number experiment. If it is a matrix the rows correspond
    to the markers and the columns to the samples.}
  \item{chrom}{the chromosomes (or other group identifier) from which
    the markers came.  Vector of length same as the number of rows of
    genomdat.  If one wants the chromosomes to be ordered in the
    natural order, this variable should be numeric or ordered category.} 
  \item{maploc}{the locations of marker on the genome.  Vector of length
    same as the number of rows of genomdat. This has to be numeric.}
  \item{data.type}{logratio (aCGH, ROMA, etc.) or binary (LOH).}
  \item{alpha}{significance levels for the test to accept change-points.}
  \item{nperm}{number of permutations used for p-value computation.}
  \item{window.size}{size of window used to speed up computations when
    segment size is too large.  Default is NULL (whole segment used).}
  \item{overlap}{proportion of data that overlap for adjacent windows.}
  \item{trim}{proportion of data to be trimmed for variance calculation
    for smoothing outliers and undoing splits based on SD.}
  \item{smooth.outliers}{should single point outliers be smoothed for
    logratio data.  Default is TRUE.}
  \item{smooth.region}{number of points to consider on the left and the
    right of a point to detect it as an outlier.}
  \item{outlier.SD}{the number of SDs away from the nearest point in the
    smoothing region to call a point an outlier.}
  \item{smooth.SD}{the number of SDs from the median in the smoothing
    region where a smoothed point is positioned.}
  \item{smooth.output}{should the smoothed data be returned.}
  \item{undo.splits}{A character string specifying how change-points are
    to be undone, if at all.  Default is "none".  Other choices are
    "prune", which uses a sum of squares criterion, and "sdundo", which 
    undoes splits that are not at least this many SDs apart.}
  \item{undo.prune}{the proportional increase in sum of squares allowed
    when eliminating splits if undo.splits="prune".}
  \item{undo.SD}{the number of SDs between means to keep a split if
    undo.splits="sdundo".}  
  \item{verbose}{if TRUE the print statements to monitor the program's
    progress are run.}
  }

\value{
  a list with components:

  \item{smoothed.data}{the smoothed data used for segmentation.  Only
    returned if smooth.output=TRUE.}
  \item{output}{a data frame with six columns.  Each row of the data
    frame contains a segment for which there are six variables: the
    sample id, the chromosome number, the map position of the start of
    the segment, the map position of the end of the segment, the number
    of markers in the segment, and the average value in the segment.}
}

\details{
  This function implements the cicular binary segmentation (CBS)
  algorithm of Olshen and Venkatraman (2004).  Given a set of genomic
  data, either continuous or binary, the algorithm recursively splits
  chromosomes into either two or three subsegments based on a maximum
  t-statistic.  A reference distribution, used to decided whether or not
  to split, is estimated by permutation.  Options are given to eliminate
  splits when the means of adjacent segments are not sufficiently far
  apart.  Note that after the first split the $\alpha$-levels of the
  tests for splitting are not unconditional.

  We recommend using one of the undoing options to remove change-points
  detected due to local trends (see the manuscript below for examples of
  local trends).

  Since the segmentation procedure uses a permutation reference
  distribution, R commands for setting and saving seeds should be used
  if the user wishes to reproduce the results.
}

\examples{

# test code on an easy data set
set.seed(25)
genomdat <- rnorm(500, sd=0.1) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
plot(genomdat)
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test1 <- segment(genomdat, chrom, maploc)

# test code on a noisier and hence more difficult data set
set.seed(51)
genomdat <- rnorm(500, sd=0.2) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
plot(genomdat)
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test2 <- segment(genomdat, chrom, maploc)

}

\author{ E. S. Venkatraman and Adam Olshen \email{olshena@mskcc.org} }

\references{
  Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004).
  Circular binary segmentation for the analysis of array-based DNA copy
  number data.  To appear in \emph{Biostatistics}.
  \url{http://www.mskcc.org/biostat/~olshena/research.}  
}


\keyword{nonparametric}
