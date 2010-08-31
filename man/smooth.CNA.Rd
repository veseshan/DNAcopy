\name{smooth.CNA}
\alias{smooth.CNA}
\alias{unsmooth.CNA}
\title{Smooth a `Copy Number Array' data object}
\description{
  Detect outliers and smooth the data prior to analysis by programs such
  as circular binary segmentation (CBS). 
}
\usage{
  smooth.CNA(x, smooth.region=10, outlier.SD.scale=4, smooth.SD.scale=2,
                trim=0.025, overwrite = FALSE)
  unsmooth.CNA(x, genomdatdiff) 
}                          
\arguments{
  \item{x}{Copy number array data object}
  \item{smooth.region}{number of points to consider on the left and the
    right of a point to detect it as an outlier.}
  \item{outlier.SD.scale}{the number of SDs away from the nearest point
    in the smoothing region to call a point an outlier.}
  \item{smooth.SD.scale}{the number of SDs from the median in the
    smoothing region where a smoothed point is positioned.}
  \item{trim}{proportion of data to be trimmed for variance calculation
    for smoothing outliers and undoing splits based on SD.}
  \item{overwrite}{when the CNA object uses file based genomdat this
    option is used to overwrite the original file.}
  \item{genomdatdiff}{a list with a matrix for each sample containing
    the locations of the probes that were modified and the original
    value. This is the output of smooth.CNA for big.matrix genomdat.}
}
\value{If the genomdat in the \code{CNA} object is a matrix then an
  object of class \code{CNA} with outliers smoothed is returned. If
  genomdat is a big.matrix then the smoothing is done in place and a
  list of locations of differences and original value.
}
\examples{

data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

#We generally recommend smoothing single point outliers before analysis
#Make sure to check that the smoothing is proper

smoothed.CNA.object <- smooth.CNA(CNA.object)

}

\keyword{nonparametric}
