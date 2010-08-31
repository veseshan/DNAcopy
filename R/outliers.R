# outlier smoothing for CNA object
smooth.CNA <- function(x, smooth.region=10, outlier.SD.scale=4,
                       smooth.SD.scale=2, trim=0.025, overwrite=FALSE) {
  if (class(x) != "CNA") stop("First arg must be of class CNA")
  if (x@data.type == "binary") stop("binary data doesn't need smoothing")

  if (!is(x@genomdat, "matrix")) {
    if (!missing(overwrite) && is.logical(overwrite) && overwrite) {
      warning("backing file of genomdat may be overwritten")
    } else {
      stop("genomdat is file backed; use overwrite=TRUE")
    }
  }

  if (overwrite) genomdatdiff <- list()
  nsample <- ncol(x@genomdat)
  chrom <- x$chrom
  uchrom <- unique(chrom)
  for (isamp in 1:nsample) {
    genomdat <- x[,isamp]
    ina <- which(is.finite(genomdat))
    trimmed.SD <- sqrt(trimmed.variance(genomdat[ina], trim))
    outlier.SD <- outlier.SD.scale*trimmed.SD
    smooth.SD <- smooth.SD.scale*trimmed.SD
    k <- smooth.region
    for (i in uchrom) {
      ina <- which(is.finite(genomdat) & chrom==i)
      n <- length(genomdat[ina])
      smoothed.data <- .Fortran("smoothLR",
                                as.integer(n),
                                as.double(genomdat[ina]),
                                sgdat=double(n),
                                as.integer(k),
                                as.double(outlier.SD),
                                as.double(smooth.SD),
                                PACKAGE = "DNAcopy")$sgdat
      x[,isamp][ina] <- smoothed.data
    }
    if (overwrite) {
      ina <- which(is.finite(genomdat))
      ii <- ina[which(genomdat[ina] != x@genomdat[ina,isamp])]
      if (length(ii) > 0) genomdatdiff[[paste("sample",isamp)]] <- cbind(loc=ii, orig=genomdat[ii])
    }
  }

  if (overwrite) {
    return(genomdatdiff)
  } else {
    return(x)
  }
}

trimmed.variance <- function(genomdat, trim=0.025) {
  n <- length(genomdat)
  n.keep <- round((1-2*trim)*(n-1))
  inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
}

inflfact <- function(trim) {
  a <- qnorm(1-trim)
  x <- seq(-a,a,length=10001)
  x1 <- (x[-10001] + x[-1])/2
  1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
}

# restore the original data when using big.matrix
unsmooth.CNA <- function(x, genomdatdiff) {
  if (class(x) != "CNA") stop("First arg must be of class CNA")
  if (!is(x@genomdat, "matrix") && length(genomdatdiff) > 0) {
    message("original logratio data being restored in backingfile")
    nn <- length(names(genomdatdiff))
    for (i in 1:nn) {
      tmp <- strsplit(names(genomdatdiff)[i], " ")[[1]]
      isamp <- as.numeric(tmp[2])
      loc <- genomdatdiff[[i]][,1]
      orig <- genomdatdiff[[i]][,2]
      x@genomdat[loc,isamp] <- orig
    }
  } else {
    message("nothing to restore")
  }
}
