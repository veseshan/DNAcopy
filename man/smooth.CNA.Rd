\name{smooth.CNA}
\alias{smooth.CNA}
\title{Smooth a `Copy Number Array' data object}
\description{
  Detect outliers and smooth the data prior to analysis by programs such
  as circular binary segmentation (CBS). 
}
\usage{
  smooth.CNA(x, smooth.region=10, outlier.SD.scale=4, smooth.SD.scale=2,
                trim=0.025)
}                          
\arguments{
  \item{x}{Copy number array data object}
  \item{smooth.region}{number of points to consider on the left and the
    right of a point to detect it as an outlier. (default=10)}
  \item{outlier.SD.scale}{the number of SDs away from the nearest point
    in the smoothing region to call a point an outlier.}
  \item{smooth.SD.scale}{the number of SDs from the median in the
    smoothing region where a smoothed point is positioned.}
  \item{trim}{proportion of data to be trimmed for variance calculation
    for smoothing outliers and undoing splits based on SD.}
}
\value{
  An object of class \code{CNA} with outliers smoothed i.e the logratio
  values of singleton outliers is shrunk towards the values of its
  neighbors. The output is of the same dimension as the input.
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
