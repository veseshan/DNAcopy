changepoints <- function(genomdat, data.type="logratio", alpha=0.01, sbdry,
                         sbn, nperm=10000, p.method="hybrid", window.size=NULL,
                         overlap=0.25, kmax=25, nmin=200, trimmed.SD = NULL,
                         undo.splits="none", undo.prune=0.05, undo.SD=3,
                         verbose=1, ngrid=100, tol=1e-6)
  {
    n <- length(genomdat)
    if (missing(trimmed.SD)) trimmed.SD <- mad(diff(genomdat))/sqrt(2)
    if (is.null(window.size)) window.size <- n
    seg.end <- c(0,n)
    k <- length(seg.end)
    change.loc <- NULL
    while (k > 1)
      {
        current.n <- seg.end[k]-seg.end[k-1]
        if (verbose>=3) cat(".... current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
        if(current.n >= 4)
          {
            current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k]]
            wsize <- min(current.n, window.size)
            winnum <- ceiling((current.n - wsize)/((1-overlap)*wsize)) + 1
            winloc <- round(seq(0,current.n - wsize,length=winnum))
            hybrid <- FALSE
            delta <- 0
            if ((p.method=="hybrid") & (nmin < current.n)) {
              hybrid <- TRUE
              delta <- (kmax+1)/current.n
            }
            zzz <- .Fortran("fndcpt",
                            n=as.integer(current.n),
                            w=as.integer(wsize),
                            twon=as.integer(2*wsize),
                            wn=as.integer(winnum),
                            wloc=as.integer(winloc),
                            x=as.double(current.genomdat),
                            px=double(current.n),
                            sx=double(wsize),
                            tx=double(2*wsize),
                            nperm=as.integer(nperm),
                            cpval=as.double(alpha),
                            ncpt=integer(1),
                            icpt=integer(2),
                            ibin=as.logical(data.type=="binary"),
                            hybrid=as.logical(hybrid),
                            hk=as.integer(kmax),
                            delta=as.double(delta),
                            ngrid=as.integer(ngrid),
                            sbn=as.integer(sbn),
                            sbdry=as.integer(sbdry),
                            tol= as.double(tol),
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
        if(verbose>=3) cat(".... segments to go:",seg.end,"\n")
      }
    seg.ends <- rev(change.loc)
    nseg <- length(seg.ends)
    lseg <- diff(c(0,seg.ends))
    if (nseg > 1) {
        if (undo.splits == "prune") {
            lseg <- changepoints.prune(genomdat, lseg, undo.prune)
        }
        if (undo.splits == "sdundo") {
            lseg <- changepoints.sdundo(genomdat, lseg, trimmed.SD, undo.SD)
        }
    }
    segmeans <- 0*lseg
    ll <- uu <- 0
    for(i in 1:length(lseg)) {
      uu <- uu + lseg[i]
      segmeans[i] <- mean(genomdat[(ll+1):uu])
      ll <- uu
    }
    list("lseg" = lseg, "segmeans" = segmeans)
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

changepoints.sdundo <- function(genomdat, lseg, trimmed.SD, change.SD=3) {
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
