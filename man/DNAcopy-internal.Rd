\name{DNAcopy-internal}
\alias{changepoints}
\alias{changepoints.prune}
\alias{changepoints.sdundo}
\alias{smooth.data}
\alias{trimmed.variance}
\alias{inflfact}
\title{Internal DNAcopy functions}
\description{
  Internal functions of package DNAcopy.
}
\usage{
changepoints(genomdat, data.type = "logratio", alpha = 0.01, nperm
                 = 10000, window.size = NULL, overlap = 0.25, 
                 trimmed.SD = NULL, undo.splits = "none", 
                 undo.prune = 0.05, undo.SD = 3, verbose = 1)
changepoints.prune(genomdat, lseg, change.cutoff=0.05)
changepoints.sdundo(genomdat, lseg, trimmed.SD, change.SD=3)
trimmed.variance(genomdat, trim=0.025)
inflfact(trim)
}
\details{These are not to be called directly by the user}
\keyword{internal}


