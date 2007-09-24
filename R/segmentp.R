segments.p <- function(x, ngrid=100, tol=1e-6)
  {
    if (!inherits(x, "DNAcopy")) 
      stop("First arg must be the result of segment")
    xdat <- x$data
    xout <- x$output
    nsample <- ncol(xdat)-2
    sampleid <- colnames(xdat)[-(1:2)]
    chrom0 <- xdat$chrom
    uchrom <- unique(chrom0)
    nchrom <- length(uchrom)
    bstat <- pval <- rep(NA, nrow(xout))
    ll <- 0
    iisamp <- 2
    for (isamp in sampleid) {
      iisamp <- iisamp + 1
      genomdat <- xdat[, iisamp]
      ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
      genomdat <- genomdat[ina]
      chrom <- chrom0[ina]
      for(ichrom in uchrom) {
        seglen <- xout$num.mark[xout$ID == isamp & xout$chrom == ichrom]
        kk <- length(seglen)
        if (kk > 1) {
          lo <- 1
          hi <- sum(seglen[1:2])
          ibstat <- ipval <- rep(NA, kk)
          gendat <- genomdat[chrom == ichrom]
          segmean <- xout$seg.mean[xout$ID == isamp & xout$chrom == ichrom]
          xresid <- gendat - rep(segmean, seglen)
          for(i in 1:(kk-1)) {
            gendati <- gendat[lo:hi]
            xresidi <- xresid[lo:hi]
            gendati <- (gendati - mean(gendati))/sd(xresidi)
            n <- length(gendati)
            zzz <- .Fortran("bsegp",
                            as.integer(n),
                            as.double(gendati),
                            ostat=double(1),
                            pval=double(1),
                            as.integer(ngrid),
                            as.double(tol),
                            PACKAGE="DNAcopy")
            ibstat[i] <- zzz$ostat
            ipval[i] <- zzz$pval
            lo <- lo + seglen[i]
            hi <- hi + seglen[i+2]
          }
          ibstat[kk] <- ipval[kk] <- NA
        } else {
          ibstat <- ipval <- NA
        }
        bstat[ll + (1:kk)] <- ibstat
        pval[ll + (1:kk)] <- ipval
        ll <- ll + kk
      }
    }
    cbind(xout, bstat, pval)
  }
