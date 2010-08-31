oldDNAcopyPlot <- function (x, plot.type=c("whole", "plateau", "samplebychrom",
                                 "chrombysample"), xmaploc=FALSE, altcol=TRUE,
                            sbyc.layout=NULL, cbys.nchrom=1, cbys.layout=NULL,
                            include.means=TRUE, zeroline=TRUE, pt.pch=NULL,
                            pt.cex=NULL, pt.cols=NULL, segcol=NULL, zlcol=NULL,
                            ylim=NULL, lwd=NULL, ...) {
  if (class(x) != "DNAcopy") stop("First arg must be the result of segment")
#  xdat <- x$data
  nsample <- ncol(x@genomdat)
#  if(missing(ylim)) {
#    uylim <- max(abs(xdat[,-(1:2)]), na.rm=TRUE)
#    ylim <- c(-uylim, uylim)
#  }
  xres <- x$output
  if(dev.cur() <= 1) dev.new()
  int.dev <- dev.interactive()
  plot.type <- match.arg(plot.type)
  op <- par(no.readonly = TRUE)
  parask <- par("ask")
  if (int.dev & !parask & nsample>1) par(ask = TRUE)
  sampleid <- colnames(x@genomdat)
  chrom0 <- x$chrom
  uchrom <- unique(chrom0)
  nchrom <- length(uchrom)
  if (xmaploc) {
    maploc0 <- as.numeric(x$maploc)
    if(length(uchrom)>1 & max(maploc0[chrom0==uchrom[1]]) > min(maploc0[chrom0==uchrom[2]])) {
      plen <- max(maploc0[chrom0==uchrom[1]])
      for(i in 2:nchrom) {
        maploc0[chrom0==uchrom[i]] <- plen + maploc0[chrom0==uchrom[i]]
        plen <- max(maploc0[chrom0==uchrom[i]])
      }
    }
  }
  if (missing(pt.pch)) pt.pch <- "."
  if (missing(pt.cex)) {
    if (pt.pch==".") { pt.cex <- 3}
    else {pt.cex <- 1}
  }
  wcol0 <- rep(1, length(chrom0))
  if (altcol) {
    j <- 0
    for (i in uchrom) {
      j <- (j+1) %% 2
      wcol0[chrom0==i] <- 1+j
    }
  }
  if (missing(pt.cols)) pt.cols <- c("black","green")
  if (missing(segcol)) segcol <- "red"
  if (missing(zlcol)) zlcol <- "grey"
  if (missing(lwd)) lwd <- 3
  if (plot.type == "chrombysample") {
    cat("Setting multi-figure configuration\n")
    par(mar = c(0, 4, 0, 2), oma = c(4, 0, 4, 0), mgp = c(2, 0.7, 0))
    if (missing(cbys.layout)) {
      nrow <- ncol <- ceiling(sqrt(nsample))
      if (nrow*ncol - nsample > 0) {
        nrow <- nrow - 1
        ncol <- ncol + 1
      }
      if (nrow*ncol - nsample >= nrow) ncol <- ncol - 1
      cbys.layout <- c(nrow, ncol)
    }
    lmat0 <- lmat1 <- c(1:nsample, rep(-cbys.nchrom*nsample, prod(cbys.layout) - nsample))
    for(i in 1:(cbys.nchrom-1)) {
      lmat1 <- c(lmat1,lmat0+nsample*i)
    }
    lmat1[lmat1<0] <- 0
    lmat <- matrix(lmat1, nrow = cbys.layout[1], ncol = cbys.nchrom*cbys.layout[2], byrow = FALSE)
    layout(lmat)
  }
  if (plot.type == "samplebychrom") {
    cat("Setting multi-figure configuration\n")
    par(mar = c(4, 4, 4, 2), oma = c(0, 0, 2, 0), mgp = c(2, 0.7, 0))
    if (missing(sbyc.layout)) {
      nrow <- ncol <- ceiling(sqrt(nchrom))
      if (nrow*ncol - nchrom > 0) {
        nrow <- nrow - 1
        ncol <- ncol + 1
      }
      if (nrow*ncol - nchrom > ncol) ncol <- ncol - 1
      sbyc.layout <- c(nrow, ncol)
    }
    lmat <- matrix(c(1:nchrom, rep(0,prod(sbyc.layout)-nchrom)),
                   nrow = sbyc.layout[1], ncol = sbyc.layout[2], byrow=TRUE)
    layout(lmat)
  }
  if (plot.type == "chrombysample") {
    atchrom <- 0.5/cbys.nchrom
    for (ichrom in uchrom) {
      for (isamp in 1:nsample) {
        genomdat <- x[chrom0==ichrom, isamp]
        ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
        genomdat <- genomdat[ina]
        ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp] & xres$chrom==ichrom]))
        mm <- xres$seg.mean[xres$ID == sampleid[isamp] & xres$chrom==ichrom]
        kk <- length(ii)
        zz <- cbind(ii[-kk] + 1, ii[-1])
        plot(genomdat, pch = pt.pch, cex=pt.cex, xaxt="n", ylim = ylim, ylab = sampleid[isamp])
        if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
        if (isamp%%cbys.layout[1] == 0) {
          axis(1, outer=TRUE)
          title(xlab="Index")
        }
        if (include.means) {
          for (i in 1:(kk - 1))
            {
              lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
            }
        }
      }
      mtext(paste("Chromosome",ichrom), side = 3, line = 1, at = atchrom, outer=TRUE, font=2)
      atchrom <- atchrom + 1/cbys.nchrom
      atchrom <- atchrom - floor(atchrom)
    }
  } else {
    for (isamp in 1:nsample)
      {
        genomdat <- x[, isamp]
        ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
        genomdat <- genomdat[ina]
        wcol <- wcol0[ina]
        chrom <- chrom0[ina]
        if (xmaploc) maploc <- maploc0[ina]
        ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]]))
        mm <- xres$seg.mean[xres$ID == sampleid[isamp]]
        kk <- length(ii)
        zz <- cbind(ii[-kk] + 1, ii[-1])
        if(missing(ylim)) ylim <- range(c(genomdat, -genomdat))
        if (plot.type=="whole")
          {
            if (xmaploc) {
              plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, col=pt.cols[wcol], main = sampleid[isamp], ylab = "", ylim = ylim)
              if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
            } else {
              plot(genomdat, pch = pt.pch, cex=pt.cex, col=pt.cols[wcol], main = sampleid[isamp], ylab = "", ylim = ylim)
              if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
            }
            if (include.means) {
              for (i in 1:(kk - 1))
                {
                  if (xmaploc) { 
                    lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
                  } else {
                    lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
                  }
                }
            }
          }
        if (plot.type=="samplebychrom")
          {
            cc <- xres$chrom[xres$ID == sampleid[isamp]]
            for (ichrom in uchrom)
              {
                plot(genomdat[chrom == ichrom], pch = pt.pch, cex=pt.cex, ylab = "", main = paste("Chromosome", ichrom), ylim = ylim)
                if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
                if (include.means) {
                  jj <- which(cc==ichrom)
                  jj0 <- min(jj)
                  for (i in jj)
                    {
                      lines(1+zz[i, ]-zz[jj0,1], rep(mm[i], 2), col = segcol, lwd=lwd)
                    }
                }
              }
            mtext(sampleid[isamp], side = 3, line = 0, outer = TRUE, font=2)
          }
        if (plot.type=="plateau")
          {
            omm <- order(mm)
            ozz <- zz[omm,]
            ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
            plot(genomdat[ina], pch = pt.pch, cex=pt.cex, main = sampleid[isamp], ylab = "", ylim = ylim)
            if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
            if (include.means) {
              ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]][omm]))
              smm <- mm[omm]
              zz <- cbind(ii[-kk] + 1, ii[-1])
              for (i in 1:(kk-1)) lines(zz[i, ], rep(smm[i], 2), col = segcol, lwd=lwd)
            }
          }
      }
  }
  on.exit( if (plot.type=="chrombysample" | plot.type=="samplebychrom") {
    par(op)
  } else { if(int.dev & !parask & nsample>1) par(ask=parask) })
}

