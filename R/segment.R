segment <- function(x, weights=NULL, alpha=0.01, nperm=10000, p.method= 
                    c("hybrid","perm"), min.width=2, kmax=25, nmin=200, 
                    eta=0.05, sbdry=NULL, trim = 0.025, undo.splits=
                    c("none","prune", "sdundo"), undo.prune=0.05, undo.SD=3,
                    verbose=1)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be a copy number array object")
    call <- match.call()
    if (min.width < 2 | min.width > 5) stop("minimum segment width should be between 2 and 5")
    if (nmin < 4*kmax) stop("nmin should be >= 4*kmax")
    if (missing(sbdry)) {
      if (nperm==10000 & alpha==0.01 & eta==0.05) {
        if (!exists("default.DNAcopy.bdry")) data(default.DNAcopy.bdry, package="DNAcopy",envir=environment())
        sbdry <- get("default.DNAcopy.bdry", envir=environment())
      } else {
        max.ones <- floor(nperm*alpha) + 1
        sbdry <- getbdry(eta, nperm, max.ones)
      }
    }
    weighted <- ifelse(missing(weights), FALSE, TRUE)
#   rudimentary error checking for weights
    if (weighted) {
      if (length(weights) != nrow(x)) stop("length of weights should be the same as the number of probes")
      if (min(weights) <= 0) stop("all weights should be positive")
    }
    sbn <- length(sbdry)
    nsample <- ncol(x)-2
    sampleid <- colnames(x)[-(1:2)]
    uchrom <- unique(x$chrom)
    data.type <- attr(x, "data.type")
    p.method <- match.arg(p.method)
    undo.splits <- match.arg(undo.splits)
    segres <- list()
    segres$data <- x
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    segRows <- list()
    segRows$startRow <- NULL
    segRows$endRow <- NULL
    for (isamp in 1:nsample) {
      if (verbose>=1) cat(paste("Analyzing:", sampleid[isamp],"\n"))
      genomdati <- x[,isamp+2]
      ina <- which(is.finite(genomdati))
      genomdati <- genomdati[ina]
      trimmed.SD <- sqrt(trimmed.variance(genomdati, trim))
      chromi <- x$chrom[ina]
#      maploci <- x$maploc[ina]
      if (weighted) {
        wghts <- weights[ina]
      } else {
        wghts <- NULL
      }
      sample.lsegs <- NULL
      sample.segmeans <- NULL
      for (ic in uchrom) {
        if (verbose>=2) cat(paste("  current chromosome:", ic, "\n"))
        segci <- changepoints(genomdati[chromi==ic], data.type, alpha, wghts,
                              sbdry, sbn, nperm, p.method, min.width, kmax,
                              nmin, trimmed.SD, undo.splits, undo.prune,
                              undo.SD, verbose)
        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)
      }
      sample.nseg <- length(sample.lsegs)
      sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
      sample.segs.end <- ina[cumsum(sample.lsegs)]
      allsegs$ID <- c(allsegs$ID, rep(isamp,sample.nseg))
      allsegs$chrom <- c(allsegs$chrom, x$chrom[sample.segs.end])
      allsegs$loc.start <- c(allsegs$loc.start, x$maploc[sample.segs.start])
      allsegs$loc.end <- c(allsegs$loc.end, x$maploc[sample.segs.end])
      allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
      allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
      segRows$startRow <- c(segRows$startRow, sample.segs.start)
      segRows$endRow <- c(segRows$endRow, sample.segs.end)
    }
    allsegs$ID <- sampleid[allsegs$ID]
    allsegs$seg.mean <- round(allsegs$seg.mean, 4)
    allsegs <- as.data.frame(allsegs)
    allsegs$ID <- as.character(allsegs$ID)
    segres$output <- allsegs
    segres$segRows <- as.data.frame(segRows)
    segres$call <- call
    if (weighted) segres$weights <- weights
    class(segres) <- "DNAcopy"
    segres
  }
