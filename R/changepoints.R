changepoints <- function(genomdat, data.type="logratio", alpha=0.01,
                         nperm=10000, window.size=NULL, overlap=0.25,
                         trim = 0.025, smooth.outliers=TRUE, smooth.region=2,
                         outlier.SD=4, smooth.SD=2, smooth.output=FALSE,
                         undo.splits="none", undo.prune=0.05, undo.SD=3,
                         verbose=TRUE)
  {
    n <- length(genomdat)
    genomdat.orig <- genomdat
    if(smooth.outliers) genomdat <- smooth.data(genomdat.orig, smooth.region, outlier.SD, smooth.SD, trim)
    if(is.null(window.size)) window.size <- n
    seg.end <- c(0,n)
    k <- length(seg.end)
    change.loc <- NULL
    while (k > 1)
      {
        current.n <- seg.end[k]-seg.end[k-1]
        if(verbose) cat("Current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
        if(current.n >= 4)
          {
            current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k]]
            wsize <- min(current.n, window.size)
            winnum <- ceiling((current.n - wsize)/((1-overlap)*wsize)) + 1
            winloc <- round(seq(0,current.n - wsize,length=winnum))
            zzz <- .Fortran("fndcpt",
                            n=as.integer(current.n),
                            w=as.integer(wsize),
                            wn=as.integer(winnum),
                            wloc=as.integer(winloc),
                            x=as.double(current.genomdat),
                            px=double(current.n),
                            sx=double(wsize),
                            tx=double(wsize),
                            nperm=as.integer(nperm),
                            cpval=as.double(alpha),
                            ncpt=integer(1),
                            icpt=integer(2),
                            ibin=as.logical(data.type=="binary"),
                            PACKAGE="DNAcopy")
          }
        else
          {
            zzz <- list()
            zzz$ncpt <- 0
          }
        if(zzz$ncpt==0) change.loc <- c(change.loc,seg.end[k])
        seg.end <- switch(1+zzz$ncpt,seg.end[-k],
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt[1],seg.end[k]),
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt,seg.end[k]))
        k <- length(seg.end)
        if(verbose) cat("Segments to go:",seg.end,"\n")
      }
    seg.ends <- rev(change.loc)
    nseg <- length(seg.ends)
    lseg <- diff(c(0,seg.ends))
    if (nseg > 1) {
        if (undo.splits == "prune") {
            lseg <- changepoints.prune(genomdat, lseg, undo.prune)
        }
        if (undo.splits == "sdundo") {
            lseg <- changepoints.sdundo(genomdat, lseg, change.SD, trim)
        }
    }
    segmeans <- 0*lseg
    if (!smooth.output) genomdat <- genomdat.orig
    ll <- uu <- 0
    for(i in 1:length(lseg)) {
      uu <- uu + lseg[i]
      segmeans[i] <- mean(genomdat[(ll+1):uu])
      ll <- uu
    }
    if (smooth.output) {
      list("smoothed.data"=genomdat, "lseg" = lseg, "segmeans" = segmeans)
    } else {
      list("lseg" = lseg, "segmeans" = segmeans)
    }
  }

changepoints.prune <- function(genomdat, lseg, change.cutoff=0.05) {
  n <- length(genomdat)
  nseg <- length(lseg)
  ncpt <- nseg-1
  zzz <- .Fortran("prune",
                  as.integer(n),
                  as.double(genomdat),
                  as.integer(nseg),
                  as.integer(lseg),
                  as.double(change.cutoff),
                  double(nseg),
                  as.integer(ncpt),
                  loc=integer(ncpt),
                  integer(2*ncpt),
                  pncpt=integer(1), PACKAGE="DNAcopy")
  pruned.ncpt <- zzz$pncpt
  pruned.cpts <- cumsum(lseg)[zzz$loc[1:pruned.ncpt]]
  pruned.lseg <- diff(c(0,pruned.cpts,n))
  pruned.lseg
}

changepoints.sdundo <- function(genomdat, lseg, change.SD=3, trim=0.025) {
  trimmed.SD <- sqrt(trimmed.variance(genomdat, trim))
  change.SD <- trimmed.SD*change.SD
  cpt.loc <- cumsum(lseg)
  sdundo <- TRUE
  while(sdundo) {
    k <- length(cpt.loc)
    if (k>1) {
      segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
      segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]])}, genomdat)
      adsegmed <- abs(diff(segmed))
      if (min(adsegmed) < change.SD) {
        i <- which(adsegmed == min(adsegmed))
        cpt.loc <- cpt.loc[-i]
      } else {
        sdundo <- FALSE
      }
    } else {
      sdundo <- FALSE
    }
  }
  lseg.sdundo <- diff(c(0,cpt.loc))
  lseg.sdundo
}

trimmed.variance <- function(genomdat, trim=0.025)
  {
    n <- length(genomdat)
    n.keep <- round((1-2*trim)*(n-1))
    inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
  }

inflfact <- function(trim)
  {
    a <- qnorm(1-trim)
    x <- seq(-a,a,length=10001)
    x1 <- (x[-10001] + x[-1])/2
    1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
  }
smooth.data <- function(genomdat, smooth.region=2, outlier.SD=4, smooth.SD=2,
                        trim=0.025)
  {
    trimmed.SD <- sqrt(trimmed.variance(genomdat, trim))
    outlier.SD <- outlier.SD*trimmed.SD
    smooth.SD <- smooth.SD*trimmed.SD
    k <- smooth.region
    n <- length(genomdat)
    smoothed.data <- sapply(1:n, function(i, x, n, nbhd, oSD, sSD) {
      xi <- x[i]
      nbhd <- i+nbhd
      xnbhd <- x[nbhd[nbhd>0 & nbhd <=n]]
      if (xi > max(xnbhd) + oSD) xi <- median(c(xi,xnbhd)) + sSD
      if (xi < min(xnbhd) - oSD) xi <- median(c(xi,xnbhd)) - sSD
      xi
    }, genomdat, n, c(-k:-1, 1:k), outlier.SD, smooth.SD)
    smoothed.data
  }

