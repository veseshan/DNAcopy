\name{convertDNAcopy}
\alias{convertCNA}
\alias{convertDNAcopy}
\title{Convert `Copy Number Array' data and results to S4 classes}
\description{
  Creates a CNA and DNAcopy data objects to new S4 classes.
}
\usage{
  convertDNAcopy(object, fbdata=NULL, backingfile=NULL, fbtype=c("ff","bm"))
  convertCNA(object, fbdata=NULL, backingfile=NULL, fbtype=c("ff","bm"))
}

\arguments{
  \item{object}{the old CNA or DNAcopy object.}
  \item{fbdata}{logical specify whether genomdat should be a file backed
  big.matrix.}
  \item{backingfile}{the backing file to be used if fbdata is TRUE.}
  \item{fbtype}{whether ff or bigmemory (bm) is used for backing file.}
}

\value{CNA or DNAcopy objects corresponding to the new S4 classes.}
  
\details{
  If fbdata option is specified the backing file will be created in the
  working directory.  A descriptor file with extension .desc will also
  be created.
}
