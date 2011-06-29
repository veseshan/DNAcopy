segments.summary <- function(x)
  {
    if (!inherits(x, "DNAcopy")) 
      stop("First arg must be the result of segment")
    xdat <- x$data
    xout <- x$output
    nsample <- ncol(xdat)-2
    sampleid <- colnames(xdat)[-(1:2)]
    seg.median <- seg.sd <- seg.mad <- rep(NA, nrow(xout))
    ll <- 0
    iisamp <- 2
    for (isamp in sampleid) {
      iisamp <- iisamp + 1
# genomdat = logratio data of sample isamp
      genomdat <- xdat[, iisamp]
# ina = location of the missing values and infinity
      ina <- which(is.finite(genomdat))
# subset out the missing & infinity locations
      genomdat <- genomdat[ina]
      seglen <- xout$num.mark[xout$ID == isamp]
      kk <- length(seglen)
      seg.sd[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), sd)
      seg.median[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), median)
      seg.mad[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), mad)
      ll <- ll + kk
    }
    xout$seg.sd <- round(seg.sd, 4)
    xout$seg.median <- round(seg.median, 4)
    xout$seg.mad <- round(seg.mad, 4)
    xout
  }
