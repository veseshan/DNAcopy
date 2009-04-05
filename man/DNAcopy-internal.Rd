\name{DNAcopy-internal}
\alias{changepoints}
\alias{changepoints.prune}
\alias{changepoints.sdundo}
\alias{smooth.data}
\alias{trimmed.variance}
\alias{inflfact}
\alias{exon.changepoint}
\title{Internal DNAcopy functions}
\description{
  Internal functions of package DNAcopy.
}
\usage{
changepoints(genomdat, data.type = "logratio", alpha = 0.01, weights = NULL,
                 sbdry, sbn, nperm = 10000, p.method="hybrid", min.width=2,
                 kmax=25, nmin = 200, trimmed.SD = NULL, undo.splits = "none",
                 undo.prune = 0.05, undo.SD = 3, verbose = 1, ngrid=100,
                 tol=1e-6) 
changepoints.prune(genomdat, lseg, change.cutoff=0.05)
changepoints.sdundo(genomdat, lseg, trimmed.SD, change.SD=3)
trimmed.variance(genomdat, trim=0.025)
inflfact(trim)
exon.changepoint(exondat, ngrid=100, tol=1e-6)
}
\details{These are not to be called directly by the user}
\keyword{internal}
