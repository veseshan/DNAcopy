\name{smooth.CNA}
\alias{smooth.CNA}
\title{Smooth a `Copy Number Array' data object}
\description{
  Detect outliers and smooth the data prior to analysis by programs such
  as circular binary segmentation (CBS). 
}
\usage{
  smooth.CNA(x, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2,
                trim=0.025)
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
}
\value{
  An object of class \code{CNA} with outliers smoothed
}
\keyword{nonparametric}
