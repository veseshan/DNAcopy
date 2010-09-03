# function to convert old CNA object
convertCNA <- function(object, fbdata=NULL, backingfile=NULL, fbtype=c("ff","bm")) {
  if (isS4(object)) {message("Already an S4 CNA object")}
  else {
    fbtype <- match.arg(fbtype)
    ns <- ncol(object) - 2
    nr <- nrow(object)
    sampleid <- colnames(object)[-(1:2)]
    chrom <- object$chrom
    if (is.numeric(chrom)) chrom <- unclass(chrom)
    if (is.factor(chrom) | is.character(chrom)) chrom <- orderedChrom(chrom)
    maploc <- object$maploc
    data.type <- attr(object, "data.type")
    if (!missing(fbdata) && is.logical(fbdata) && fbdata) {
      if (fbtype == "ff") {
        genomdat <- ff(vmode="double", dim=c(nr, ns), file=backingfile)
        colnames(genomdat) <- sampleid
      } else {
        genomdat <- filebacked.big.matrix(nr, ns, backingfile=backingfile, dimnames=list(c(),sampleid))
      }
    } else {
      fbdata <- FALSE
      genomdat <- matrix(0, nr,ns)
      colnames(genomdat) <- sampleid
    }
    for (i in 1:ns) genomdat[,i] <- object[,i+2]
    if (fbtype=="bm") {
      bigmemdesc <- describe(genomdat)@description
      new("CNA", chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type, bigmemorydesc=bigmemdesc)
    } else {
      new("CNA", chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type)
    }
  }
}

# function to convert old DNAcopy object
convertDNAcopy <- function(object, fbdata=NULL, backingfile=NULL, fbtype=c("ff","bm")) {
  if (isS4(object)) {message("Already an S4 CNA object")}
  else {
    fbtype <- match.arg(fbtype)
    call <- object$call
    ns <- ncol(object$data) - 2
    nr <- nrow(object$data)
    sampleid <- colnames(object$data)[-(1:2)]
    chrom <- object$data$chrom
    if (is.numeric(chrom)) chrom <- unclass(chrom)
    if (is.factor(chrom) | is.character(chrom)) chrom <- orderedChrom(chrom)
    maploc <- object$data$maploc
    output <- object$output
    if (is.null(object$segRows)) {segRows <- data.frame()}
    else {segRows <- object$segRows}
    data.type <- attr(object$data, "data.type")
    if (!missing(fbdata) && is.logical(fbdata) && fbdata) {
      if (fbtype == "ff") {
        genomdat <- ff(vmode="double", dim=c(nr, ns), file=backingfile)
        colnames(genomdat) <- sampleid
      } else {
        genomdat <- filebacked.big.matrix(nr, ns, backingfile=backingfile, dimnames=list(c(),sampleid))
      }
    } else {
      fbdata <- FALSE
      genomdat <- matrix(0, nr,ns)
      colnames(genomdat) <- sampleid
    }
    for (i in 1:ns) genomdat[,i] <- object$data[,i+2]
    if (fbtype=="bm") {
      bigmemdesc <- describe(genomdat)@description
      new("DNAcopy", call=call, output=output, segRows= segRows, chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type, fbspecs=fbspecs)
    } else {
      new("DNAcopy", call=call, output=output, segRows= segRows, chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type)
    }
  }
}
