segment <- function(genomdat, chrom, maploc, data.type=c("logratio","binary"),
                    alpha=0.01, nperm=10000, window.size=NULL, overlap=0.25,
                    trim = 0.025, smooth.outliers=TRUE, smooth.region=2,
                    outlier.SD=4, smooth.SD=2, smooth.output=FALSE,
                    undo.splits= c("none","prune","sdundo"), undo.prune=0.05,
                    undo.SD=3, verbose=TRUE)
  {
    sortindex <- order(chrom, maploc)
    if (is.matrix(genomdat)) {
        genomdat <- genomdat[sortindex, ]
    }
    else {
        genomdat <- genomdat[sortindex]
    }
    if (is.vector(genomdat)) genomdat <- as.matrix(genomdat)
    if (smooth.output) genomdat.smoothed <- matrix(NA, nrow(genomdat), ncol(genomdat))
    nsample <- ncol(genomdat)
    chrom <- chrom[sortindex]
    maploc <- maploc[sortindex]
    uchrom <- unique(chrom)
    data.type <- match.arg(data.type)
    if(data.type=="binary") smooth.outliers <- FALSE
    undo.splits <- match.arg(undo.splits)
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    for (isamp in 1:nsample) {
      genomdati <- genomdat[,isamp]
      ina <- which(!is.na(genomdati) & !(abs(genomdati)==Inf))
      genomdati <- genomdati[ina]
      chromi <- chrom[ina]
      sample.lsegs <- NULL
      sample.segmeans <- NULL
      for (ic in uchrom) {
        if(verbose) cat(paste("Sample:", isamp, "; chrom:", ic, "\n"))
        segci <- changepoints(genomdati[chromi==ic], data.type, alpha, 
                              nperm,  window.size, overlap, trim, 
                              smooth.outliers, smooth.region, outlier.SD, 
                              smooth.SD, smooth.output, undo.splits, 
                              undo.prune, undo.SD, verbose)
        if (smooth.output) {
          genomdat.smoothed[ina, isamp][chromi==ic] <- segci$smoothed.data
        }
        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)
      }
      sample.nseg <- length(sample.lsegs)
      sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
      sample.segs.end <- ina[cumsum(sample.lsegs)]
      allsegs$ID <- c(allsegs$ID, rep(isamp,sample.nseg))
      allsegs$chrom <- c(allsegs$chrom, chrom[sample.segs.end])
      allsegs$loc.start <- c(allsegs$loc.start, maploc[sample.segs.start])
      allsegs$loc.end <- c(allsegs$loc.end, maploc[sample.segs.end])
      allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
      allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
    }
    allsegs$seg.mean <- round(allsegs$seg.mean, 4)
    if (smooth.output) {
      list(smoothed.data=genomdat.smoothed, output=as.data.frame(allsegs))
    } else {
      list(output=as.data.frame(allsegs))
    }
  }
