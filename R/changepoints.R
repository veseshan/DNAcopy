changepoints <- function(genomdat, data.type="logratio", alpha=0.01, weights=
                         NULL, sbdry, sbn, nperm=10000, p.method="hybrid", 
                         min.width=2, kmax=25, nmin=200, trimmed.SD=NULL, 
                         undo.splits="none", undo.prune=0.05, undo.SD=3,
                         verbose=1, ngrid=100, tol=1e-6)
  {
    n <- length(genomdat)
    if (missing(trimmed.SD)) trimmed.SD <- mad(diff(genomdat))/sqrt(2)
#   start with the whole 
    seg.end <- c(0,n)
    k <- length(seg.end)
    change.loc <- NULL
    weighted <- ifelse(is.null(weights), FALSE, TRUE) 
    while (k > 1)
      {
        current.n <- seg.end[k]-seg.end[k-1]
        if (verbose>=3) cat(".... current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
        if(current.n >= 2*min.width) {
          current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k]]
#   check whether hybrid method needs to be used
          hybrid <- FALSE
          delta <- 0
          if ((p.method=="hybrid") & (nmin < current.n)) {
            hybrid <- TRUE
            delta <- (kmax+1)/current.n
          }
#   call the changepoint routine
          if (weighted) {
#   get the weights for the current set of probes
            current.wts <- weights[(seg.end[k-1]+1):seg.end[k]]
            current.rwts <- sqrt(current.wts)
            current.cwts <- cumsum(current.wts)/sqrt(sum(current.wts))
#   if all values of current.genomdat are the same don't segment
            if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
              zzz <- list()
              zzz$ncpt <- 0
            } else {
#   centering the current data will save a lot of computations later
              current.avg <- sum(current.genomdat*current.wts)/sum(current.wts)
              current.genomdat <- current.genomdat - current.avg
#   need total sum of squares too
              current.tss <- sum(current.wts*(current.genomdat^2))
              zzz <- .Fortran("wfindcpt",
                              n=as.integer(current.n),
                              x=as.double(current.genomdat),
                              tss=as.double(current.tss),
                              wts=as.double(current.wts),
                              rwts=as.double(current.rwts),
                              cwts=as.double(current.cwts),
                              px=double(current.n),
                              sx=double(current.n),
                              nperm=as.integer(nperm),
                              cpval=as.double(alpha),
                              ncpt=integer(1),
                              icpt=integer(2),
                              hybrid=as.logical(hybrid),
                              al0=as.integer(min.width),
                              hk=as.integer(kmax),
                              mncwt=double(kmax),
                              delta=as.double(delta),
                              ngrid=as.integer(ngrid),
                              sbn=as.integer(sbn),
                              sbdry=as.integer(sbdry),
                              tol= as.double(tol),
                              PACKAGE="DNAcopy")
            }
          } else { 
#   if all values of current.genomdat are the same don't segment
            if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
              zzz <- list()
              zzz$ncpt <- 0
            } else {
#   centering the current data will save a lot of computations later
              current.avg <- mean(current.genomdat)
              current.genomdat <- current.genomdat - current.avg
#   need total sum of squares too
              current.tss <- sum(current.genomdat^2)
              zzz <- .Fortran("fndcpt",
                              n=as.integer(current.n),
                              x=as.double(current.genomdat),
                              tss=as.double(current.tss),
                              px=double(current.n),
                              sx=double(current.n),
                              nperm=as.integer(nperm),
                              cpval=as.double(alpha),
                              ncpt=integer(1),
                              icpt=integer(2),
                              ibin=as.logical(data.type=="binary"),
                              hybrid=as.logical(hybrid),
                              al0=as.integer(min.width),
                              hk=as.integer(kmax),
                              delta=as.double(delta),
                              ngrid=as.integer(ngrid),
                              sbn=as.integer(sbn),
                              sbdry=as.integer(sbdry),
                              tol= as.double(tol),
                              PACKAGE="DNAcopy")
            }
          }
        } else {
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
      if (weighted) {
        segmeans[i] <- sum(genomdat[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
      } else {
        segmeans[i] <- mean(genomdat[(ll+1):uu])
      }
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
    x <- seq(-a,a,length.out=10001)
    x1 <- (x[-10001] + x[-1])/2
    1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
  }
