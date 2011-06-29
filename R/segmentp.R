segments.p <- function(x, ngrid=100, tol=1e-6, alpha=0.05, search.range=100,
                       nperm=1000)
  {
    if (!inherits(x, "DNAcopy")) 
      stop("First arg must be the result of segment")
    xdat <- x$data
    xout <- x$output
    nsample <- ncol(xdat)-2
    sampleid <- colnames(xdat)[-(1:2)]
    chrom0 <- xdat$chrom
    maploc0 <- xdat$maploc
    uchrom <- unique(chrom0)
    nchrom <- length(uchrom)
    bstat <- pval <- lcl <- ucl <- rep(NA, nrow(xout))
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
      chrom <- chrom0[ina]
      maploc <- maploc0[ina]
      for(ichrom in uchrom) {
# kk = number of segments in chromosome ichrom of sample isamp
        kk <- sum(1*(xout$ID == isamp & xout$chrom == ichrom))
        if (kk > 1) {
# gendat = logratio data in chromosome ichrom of sample isamp
          gendat <- genomdat[chrom == ichrom]
# seglen = lengths of the segments in chromosome ichrom of sample isamp
          seglen <- xout$num.mark[xout$ID == isamp & xout$chrom == ichrom]
# segmean = means of the segments in chromosome ichrom of sample isamp
          segmean <- xout$seg.mean[xout$ID == isamp & xout$chrom == ichrom]
# xresid = residuals of the data in chromosome ichrom of sample isamp
          xresid <- gendat - rep(segmean, seglen)
          ibstat <- ipval <- ilcl <- iucl <- rep(NA, kk)
# begin with the first 2 segments lo & hi are the start & end points
          lo <- 1
          hi <- sum(seglen[1:2])
          for(i in 1:(kk-1)) {
# prep data from adjacent segments
            gendati <- gendat[lo:hi]
            xresidi <- xresid[lo:hi]
# standardize data
            gendati <- (gendati - mean(gendati))/sd(xresidi)
            n <- length(gendati)
# call the p-value subroutine
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
# additional data for CI routine 
# k = location of change-point
# sr = search range
# sumxk = partial sum at k (all paths are pegged at that point)
# var.factor = variance for 2-sample t-statistic
            k <- seglen[i]
            sr <- c(max(2, k-search.range),min(n-2,k+search.range))
            sumxk <- sum(gendati[1:k])
            var.factor <- n/((1:n)*(n:1 - 1))
            var.factor[n] <- 0
# call the confidence subroutine
            zzz <- .Fortran("bsegci",
                            as.integer(n),
                            as.integer(k),
                            as.double(sumxk),
                            as.double(gendati),
                            px = double(n),
                            sr = as.integer(sr),
                            vfact = as.double(var.factor),
                            as.integer(nperm),
                            bsloc = integer(nperm),
                            PACKAGE="DNAcopy")
            bsloc <- zzz$bsloc
            bsci <- quantile(bsloc, c(alpha/2, 1-alpha/2), type=1)
            ilcl[i] <- bsci[1]
            iucl[i] <- bsci[2]
# increment to the next segment
            lo <- lo + seglen[i]
            if(i < kk-1) hi <- hi + seglen[i+2]
          }
          ibstat[kk] <- ipval[kk] <- ilcl[kk] <- iucl[kk] <- NA
        } else {
          seglen <- ibstat <- ipval <- ilcl <- iucl <- NA
        }
        bstat[ll + (1:kk)] <- ibstat
        pval[ll + (1:kk)] <- ipval
# convert the lcl & ucl from probe number to maploc
        lcl[ll + (1:kk)] <- maploc[chrom == ichrom][cumsum(seglen) + (ilcl - seglen)]
        ucl[ll + (1:kk)] <- maploc[chrom == ichrom][cumsum(seglen) + (iucl - seglen)]
        ll <- ll + kk
      }
    }
    cbind(xout, bstat, pval, lcl, ucl)
  }
