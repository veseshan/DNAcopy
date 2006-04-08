\name{getbdry}
\alias{getbdry}
\alias{default.DNAcopy.bdry}
\title{Sequential stopping boundary}
\description{
   Function to compute the sequential boundary for early stopping.
}
\usage{
  getbdry(eta, nperm, max.ones, tol= 1e-2)
}                          
\arguments{
  \item{eta}{Type I error rate of the boundary.}
  \item{nperm}{Number of permutations for the reference distribution.}
  \item{max.ones}{maximum number of ones given by "floor(nperm*alpha)+1".} 
  \item{tol}{tolerance level for the iterations.}
}
\value{
  A vector integer values of length max.ones*(max.ones+1)/2 corresponding
  to the boundary for the number of ones from 1 to max.ones. The default
  boundary for nperm=10000, eta=0.05, alpha=0.01 is stored in the data
  object "default.DNAcopy.bdry".  Use this function to get the boundary
  for your favorite values for the parameters "nperm, eta, alpha" and use
  it for the argument "sbdry" in the function "segment."
}
\keyword{nonparametric}
