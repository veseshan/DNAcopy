glFrequency <- function(xout, threshold=1) {
  if (class(xout) != "DNAcopy") stop("First arg must be of class DNAcopy")
  nsample <- ncol(xout@genomdat)
  snames <- colnames(xout@genomdat)
  xmad <- rep(NA,nsample)
# initialize the probe level counts
  pfreq <- gain <- loss <- rep(0, nrow(xout@genomdat))
  for(i in 1:nsample) {
# data for ith sample
    genomdat <- xout[,i]
# ii = locations of non-missing values and finite
    ii <- which(is.finite(genomdat))
# segment means as a vector
    segout <- xout$output[xout$output$ID==snames[i],]
    segmean <- rep(segout$seg.mean, segout$num.mark)
    xmad[i] <- mad(genomdat[ii] - rep(segmean))
# gains and losses
    pfreq[ii] <- pfreq[ii] + 1
    gain[ii] <- gain[ii] + 1*((segmean - median(segmean))/xmad[i] > threshold)
    loss[ii] <- loss[ii] - 1*((segmean - median(segmean))/xmad[i] < -threshold)
  }
  out <- list()
  out$chrom <- xout$chrom
  out$maploc <- xout$maploc
  out$pfreq <- pfreq
  out$gain <- gain/pfreq
  out$loss <- loss/pfreq
  as.data.frame(out)
}
