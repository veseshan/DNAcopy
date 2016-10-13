CNA <- function(genomdat, chrom, maploc, data.type=c("logratio","binary"),
                   sampleid=NULL, presorted=FALSE)
  {
    if (is.data.frame(genomdat)) genomdat <- as.matrix(genomdat)
    if (!is.numeric(genomdat)) stop("genomdat must be numeric")
    if (!is.numeric(maploc)) stop("maploc must be numeric")
    data.type <- match.arg(data.type)
    ina <- (!is.na(chrom) & is.finite(maploc))
    if (sum(!ina)>0)
      warning("markers with missing chrom and/or maploc removed\n")
    if (!presorted) {
      sortindex <- which(ina)[order(chrom[ina], maploc[ina])]
    } else {
      sortindex <- which(ina)
    }
    if (is.factor(chrom)) chrom <- as.character(chrom)
    # added to allow arrays of single dimension - results from data.frame ops
    if (is.array(genomdat)) {
        if (length(dim(genomdat)) == 1) {
            genomdat <- as.matrix(genomdat)
        }
    }
    if (is.vector(genomdat)) genomdat <- as.matrix(genomdat)
    if (!missing(sampleid)) {
      if (length(sampleid) != ncol(genomdat)) {
        warning("length(sampleid) and ncol(genomdat) differ, names ignored\n")
        sampleid <- paste("Sample", 1:ncol(genomdat))
      } 
    } else {
        sampleid <- paste("Sample", 1:ncol(genomdat))
    }
    colnames(genomdat) <- sampleid
    zzz <- data.frame(chrom=I(chrom), maploc=maploc, genomdat)
    zzz <- zzz[sortindex,]

# check for duplicate probes (i.e. repeated maploc within a chromosome)
    if (length(ii <- which(diff(maploc)==0)) > 0) {
      if (any(chrom[ii]==chrom[ii+1])) warning("array has repeated maploc positions\n")
    }

    attr(zzz, "data.type") <- data.type
    class(zzz) <- c("CNA","data.frame")
    zzz
  }

subset.CNA <- function(x, chromlist=NULL, samplelist=NULL, ...)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be of class CNA")
    chrom <- x$chrom
    uchrom <- unique(chrom)
    if (missing(chromlist)) chromlist <- uchrom
    if (length(setdiff(chromlist, uchrom)) > 0)
      stop("chromlist contains chromosomes not in the data")
    if (length(chromlist) > length(unique(chromlist)))
      warning("duplicate chromosomes in chromlist removed")
    sampleid <- colnames(x)[-(1:2)]
    if (missing(samplelist)) samplelist <- sampleid
    nsample <- length(sampleid)
    if (length(setdiff(samplelist, 1:nsample)) > 0 & length(setdiff(samplelist, sampleid)) > 0)
      stop("samplelist should be a list of valid sample numbers or names")
    if (!is.numeric(samplelist)) samplelist <- match(samplelist, names(x)) - 2
    if (length(samplelist) > length(unique(samplelist)))
      warning("duplicate samples in samplelist removed")
    samplelist <- unique(samplelist)
    y <- x[chrom %in% chromlist,c(1:2,samplelist+2)]
    attr(y, "data.type") <- attr(x, "data.type")
    y
  }

smooth.CNA <- function(x, smooth.region=10, outlier.SD.scale=4,
                       smooth.SD.scale=2, trim=0.025)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be of class CNA")
    nsample <- ncol(x)-2
    chrom <- x$chrom
    uchrom <- unique(chrom)
    if(attr(x, "data.type")=="binary") stop("Not smoothing binary data ")
    for (isamp in 1:nsample) {
      genomdat <- x[,isamp+2]
      ina <- which(is.finite(genomdat))
      trimmed.SD <- sqrt(trimmed.variance(genomdat[ina], trim))
      outlier.SD <- outlier.SD.scale*trimmed.SD
      smooth.SD <- smooth.SD.scale*trimmed.SD
      k <- smooth.region
      n <- length(genomdat[ina])
      cfrq <- diff(c(which(!duplicated(chrom[ina])), n+1))
      nchr <- length(cfrq) # to allow for some chrom with all missing
      smoothed.data <- .Fortran("smoothLR",
                                as.integer(n),
                                as.double(genomdat[ina]),
                                as.integer(nchr),
                                as.integer(cfrq),
                                sgdat=double(n),
                                as.integer(k),
                                as.double(outlier.SD),
                                as.double(smooth.SD),
                                PACKAGE = "DNAcopy")$sgdat
      x[,isamp+2][ina] <- smoothed.data
    }
    x
  }

print.CNA <- function(x, ...)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be of class CNA")
    cat("Number of Samples", ncol(x)-2,
        "\nNumber of Probes ", nrow(x),
        "\nData Type        ", attr(x,"data.type"),"\n")
  }

plot.DNAcopy <- function (x, plot.type=c("whole", "plateau", "samplebychrom",
                               "chrombysample"), xmaploc=FALSE, altcol=TRUE,
                          sbyc.layout=NULL, cbys.nchrom=1, cbys.layout=NULL,
                          include.means=TRUE, zeroline=TRUE, pt.pch=NULL,
                          pt.cex=NULL, pt.cols=NULL, segcol=NULL, zlcol=NULL,
                          ylim=NULL, lwd=NULL, ...)
{
  if (!inherits(x, "DNAcopy")) 
    stop("First arg must be the result of segment")
  xdat <- x$data
  nsample <- ncol(xdat)-2
  if(missing(ylim)) {
    uylim <- max(abs(xdat[,-(1:2)]), na.rm=TRUE)
    ylim <- c(-uylim, uylim)
  }
  xres <- x$output
  if(dev.cur() <= 1) dev.new()
  int.dev <- dev.interactive()
  plot.type <- match.arg(plot.type)
  op <- par(no.readonly = TRUE)
  parask <- par("ask")
  if (int.dev & !parask & nsample>1) par(ask = TRUE)
  sampleid <- colnames(xdat)[-(1:2)]
  chrom0 <- xdat$chrom
  uchrom <- unique(chrom0)
  nchrom <- length(uchrom)
  if (xmaploc) {
    maploc0 <- as.numeric(xdat$maploc)
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
      if (xmaploc) maploc1 <- maploc0[chrom0==ichrom]
      for (isamp in 1:nsample) {
        genomdat <- xdat[chrom0==ichrom, isamp+2]
        ina <- which(is.finite(genomdat))
        genomdat <- genomdat[ina]
        if (xmaploc) maploc <- maploc1[ina]
        ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp] & xres$chrom==ichrom]))
        mm <- xres$seg.mean[xres$ID == sampleid[isamp] & xres$chrom==ichrom]
        kk <- length(ii)
        zz <- cbind(ii[-kk] + 1, ii[-1])
        if (xmaploc) {
          plot(maploc, genomdat, pch = pt.pch, cex=pt.cex, xaxt="n", ylim = ylim, ylab = sampleid[isamp])
        } else {
          plot(genomdat, pch = pt.pch, cex=pt.cex, xaxt="n", ylim = ylim, ylab = sampleid[isamp])
        }
        if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
        if (isamp%%cbys.layout[1] == 0) {
          axis(1, outer=TRUE)
          title(xlab="Index")
        }
        if (include.means) {
          if (xmaploc) {
            segments(maploc[zz[,1]], mm, x1=maploc[zz[,2]], y1=mm, col = segcol, lwd=lwd)
          } else {
            segments(zz[,1], mm, x1=zz[,2], y1=mm, col = segcol, lwd=lwd)
          }
#          for (i in 1:(kk - 1)) {
#            if (xmaploc) { 
#              lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
#            } else {
#              lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
#            }
#          }
        }
      }
      mtext(paste("Chromosome",ichrom), side = 3, line = 1, at = atchrom, outer=TRUE, font=2)
      atchrom <- atchrom + 1/cbys.nchrom
      atchrom <- atchrom - floor(atchrom)
    }
  } else {
    for (isamp in 1:nsample)
      {
        genomdat <- xdat[, isamp+2]
        ina <- which(is.finite(genomdat))
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
              if (xmaploc) {
                segments(maploc[zz[,1]], mm, x1=maploc[zz[,2]], y1=mm, col = segcol, lwd=lwd)
              } else {
                segments(zz[,1], mm, x1=zz[,2], y1=mm, col = segcol, lwd=lwd)
              }
#              for (i in 1:(kk - 1))
#                {
#                  if (xmaploc) { 
#                    lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
#                  } else {
#                    lines(zz[i, ], rep(mm[i], 2), col = segcol, lwd=lwd)
#                  }
#                }
            }
          }
        if (plot.type=="samplebychrom")
          {
            cc <- xres$chrom[xres$ID == sampleid[isamp]]
            for (ichrom in uchrom)
              {
                if (xmaploc) {
                  plot(maploc[chrom == ichrom], genomdat[chrom == ichrom], pch = pt.pch, cex=pt.cex, xlab="maploc", ylab = "", main = paste("Chromosome", ichrom), ylim = ylim)
                } else {
                  plot(genomdat[chrom == ichrom], pch = pt.pch, cex=pt.cex, ylab = "", main = paste("Chromosome", ichrom), ylim = ylim)
                }
                if(zeroline) abline(h=0, col=zlcol, lwd=lwd)
                if (include.means) {
                  jj <- which(cc==ichrom)
                  jj0 <- min(jj)
                  if (xmaploc) {
                    segments(maploc[zz[jj,1]], mm[jj], x1=maploc[zz[jj,2]], y1=mm[jj], col = segcol, lwd=lwd)
                  } else {
                    segments(1+zz[jj,1]-zz[jj0,1], mm[jj], x1=1+zz[jj,2]-zz[jj0,1], y1=mm[jj], col = segcol, lwd=lwd)
                  }
#                  for (i in jj)
#                    {
#                      if (xmaploc) {
#                        lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, lwd=lwd)
#                      } else {
#                        lines(1+zz[i, ]-zz[jj0,1], rep(mm[i], 2), col = segcol, lwd=lwd)
#                      }
#                    }
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
              segments(zz[,1], smm, x1=zz[,2], y1=smm, col = segcol, lwd=lwd)
#              for (i in 1:(kk-1)) lines(zz[i, ], rep(smm[i], 2), col = segcol, lwd=lwd)
            }
          }
      }
  }
  on.exit( if (plot.type=="chrombysample" | plot.type=="samplebychrom") {
    par(op)
  } else { if(int.dev & !parask & nsample>1) par(ask=parask) })
}

print.DNAcopy <- function(x, showSegRows=FALSE, ...)
  {
    if (!inherits(x, "DNAcopy")) stop("Object is not the result of segment")
    if (!is.null(cl<- x$call))
      {
        cat("Call:\n")
        dput(cl)
        cat("\n")
      }
    if (showSegRows) {
      if (is.null(x$segRows)) {
        print(x$output)
        warning("segRows missing.  Object may be a subset or from DNAcopy < 1.23.2.\n")
      } else {
        print(cbind(x$output, x$segRows))
      }
    } else {
      print(x$output)
    }
  }

subset.DNAcopy <- function(x, chromlist=NULL, samplelist=NULL, ...)
  {
    if (!inherits(x, 'DNAcopy')) stop("First arg must be of class DNAcopy")
    zdat <- x$data
    zres <- x$output
    chrom <- zdat$chrom
    uchrom <- unique(chrom)
    if (missing(chromlist) | is.null(chromlist)) chromlist <- uchrom
    if (length(setdiff(chromlist, uchrom)) > 0)
      stop("chromlist contains chromosomes not in the data")
    if (length(chromlist) > length(unique(chromlist)))
      warning("duplicate chromosomes in chromlist removed")
    sampleid <- colnames(zdat)[-(1:2)]
    if (missing(samplelist)) samplelist <- sampleid
    nsample <- length(sampleid)
    if (length(setdiff(samplelist, 1:nsample)) > 0 & length(setdiff(samplelist, sampleid)) > 0)
      stop("samplelist should be a list of valid sample numbers or names")
    if (!is.numeric(samplelist)) samplelist <- match(samplelist, names(zdat)) - 2
    if (length(samplelist) > length(unique(samplelist)))
      warning("duplicate samples in samplelist removed")
    samplelist <- unique(samplelist)
    jj <- unlist(sapply(sampleid[samplelist], function(i, id) {which(id==i)}, zres$ID ))
    zres <- zres[jj,]
    y <- list()
    y$data <- zdat[chrom %in% chromlist,c(1:2,samplelist+2)]
    attr(y$data, "data.type") <- attr(zdat, "data.type")
    y$output <- zres[zres$chrom %in% chromlist,]
    class(y) <- "DNAcopy"
    y
  }

# Chromosome.Lengths <- c(263, 255, 214, 203, 194, 183, 171, 155, 145, 144, 144, 143, 114, 109, 106, 98, 92, 85, 67, 72, 50, 56, 164, 59)
# names(Chromosome.Lengths) <- c(as.character(1:22),"X","Y")
