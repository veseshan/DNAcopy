\name{segment}
\alias{segment}
\title{Genome Segmentation Program}
\description{
  This program segments DNA copy number data into regions of estimated 
  equal copy number using circular binary segmentation (CBS).
}
\usage{
  segment(x, weights = NULL, alpha = 0.01, nperm = 10000, p.method = 
                    c("hybrid", "perm"), min.width=2, kmax=25, nmin=200, 
                    eta=0.05, sbdry=NULL, trim = 0.025, undo.splits = 
                    c("none", "prune", "sdundo"), undo.prune=0.05, 
                    undo.SD=3, verbose=1)
}
\arguments{
  \item{x}{an object of class CNA}
  \item{weights}{a vector of weights for the probes. The weights should be
    inversely proportional to their variances.  Currently all weights
    should be positive i.e. remove probes with zero weight prior to
    segmentation.}
  \item{alpha}{significance levels for the test to accept change-points.}
  \item{nperm}{number of permutations used for p-value computation.}
  \item{p.method}{method used for p-value computation.  For the "perm"
    method the p-value is based on full permutation.  For the "hybrid"
    method the maximum over the entire region is split into maximum of
    max over small segments and max over the rest.  Approximation is
    used for the larger segment max. Default is hybrid.}
  \item{min.width}{the minimum number of markers for a changed segment.
    The default is 2 but can be made larger.  Maximum possible value is
    set at 5 since arbitrary widths can have the undesirable effect of
    incorrect change-points when a true signal of narrow widths exists.}
  \item{kmax}{the maximum width of smaller segment for permutation
    in the hybrid method.}
  \item{nmin}{the minimum length of data for which the approximation of
    maximum statistic is used under the hybrid method. should be larger
    than 4*kmax}
  \item{eta}{the probability to declare a change conditioned on the
    permuted statistic exceeding the observed statistic exactly  
    j (= 1,...,nperm*alpha) times.}
  \item{sbdry}{the sequential boundary used to stop and declare a
    change. This boundary is a function of nperm, alpha and eta.  It can
    be obtained using the function "getbdry" and used instead of having
    the "segment" function compute it every time it is called.}
  \item{trim}{proportion of data to be trimmed for variance calculation
    for smoothing outliers and undoing splits based on SD.}
  \item{undo.splits}{A character string specifying how change-points are
    to be undone, if at all.  Default is "none".  Other choices are
    "prune", which uses a sum of squares criterion, and "sdundo", which 
    undoes splits that are not at least this many SDs apart.}
  \item{undo.prune}{the proportional increase in sum of squares allowed
    when eliminating splits if undo.splits="prune".}
  \item{undo.SD}{the number of SDs between means to keep a split if
    undo.splits="sdundo".}
  \item{verbose}{level of verbosity for monitoring the program's
    progress where 0 produces no printout, 1 prints the current sample,
    2 the current chromosome and 3 the current segment.  The default
    level is 1.}  
}

\value{
  An object of class \code{DNAcopy}.  It has three elements:
  \item{data}{The original CNA object which was the input for segment}
  \item{out}{a data frame with six columns.  Each row of the data frame
    contains a segment for which there are six variables: the sample id,
    the chromosome number, the map position of the start of the segment,
    the map position of the end of the segment, the number of markers in
    the segment, and the average value in the segment.}
  \item{segRows}{a data frame with the start and end row of each segment
    in the data matrix.  print command shows it with the showSegRows=T}
  \item{call}{the call that produced the output object.}
}

\details{
  This function implements the cicular binary segmentation (CBS)
  algorithm of Olshen and Venkatraman (2004).  Given a set of genomic
  data, either continuous or binary, the algorithm recursively splits
  chromosomes into either two or three subsegments based on a maximum
  t-statistic.  A reference distribution, used to decided whether or not
  to split, is estimated by permutation.  Options are given to eliminate
  splits when the means of adjacent segments are not sufficiently far
  apart.  Note that after the first split the \eqn{\alpha}-levels of the
  tests for splitting are not unconditional.

  We recommend using one of the undoing options to remove change-points
  detected due to local trends (see the manuscript below for examples of
  local trends).

  Since the segmentation procedure uses a permutation reference
  distribution, R commands for setting and saving seeds should be used
  if the user wishes to reproduce the results.

  Data that are NA, Inf, NaN will be removed on a per sample basis for
  "genomdat" and all samples for "chrom" and "maploc".  
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

# test code on a noisier and hence more difficult data set
set.seed(51)
genomdat <- rnorm(500, sd=0.2) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
plot(genomdat)
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test2 <- segment(CNA(genomdat, chrom, maploc))

# test code for weighted CBS
set.seed(97)
wts <- sample(1:3, 500, replace=TRUE)
genomdat <- rnorm(500, sd=0.3)/sqrt(wts) +
rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
plot(genomdat)
chrom <- rep(1:2,c(290,210))
maploc <- c(1:290,1:210)
test3 <- segment(CNA(genomdat, chrom, maploc), weights=wts)

#A real analyis

data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

#We generally recommend smoothing single point outliers before analysis
#Make sure to check that the smoothing is proper

smoothed.CNA.object <- smooth.CNA(CNA.object)

#Segmentation at default parameters

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

#We generally recommend smoothing single point outliers before analysis
#Make sure to check that the smoothing is proper

smoothed.CNA.object <- smooth.CNA(CNA.object)

#Segmentation at default parameters

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

}

\author{Venkatraman E. Seshan \email{seshanv@mskcc.org} and Adam Olshen
  \email{olshena@biostat.ucsf.edu} } 

\references{
  Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004).
  Circular binary segmentation for the analysis of array-based DNA copy
  number data.  \emph{Biostatistics} 5: 557-572.

  Venkatraman, E. S., Olshen, A. B. (2007) A faster circular binary
  segmentation algorithm for the analysis of array CGH data.
  \emph{Bioinformatics} 23: 657-63.
}

\keyword{nonparametric}
