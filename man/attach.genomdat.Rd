\name{attach.genomdat}
\alias{attach.genomdat}
\title{Attach genomdat in a restored CNA or DNAcopy object.}
\description{
  Creates the necessary links when filebacked genomdat is used.
}
\usage{
  attach.genomdat(obj)
}

\arguments{
  \item{obj}{the old CNA or DNAcopy object.}
}

\details{
  This function is useful when a R session where filebacked data are
  used in DNAcopy.  When the variables are restored in a future R
  session the external pointers referring to the backing file will be
  nil and thus need to be restored.  This function uses open.ff in ff
  and attach.resource in bigmemory to complete this link.

  You need to load the package ff or bigmemory prior to invoking this.
}

