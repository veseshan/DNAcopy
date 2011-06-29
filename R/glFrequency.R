glFrequency <- function(xout, threshold=1) {
  if (!inherits(xout, 'DNAcopy')) stop("First arg must be of class DNAcopy")
  nsample <- ncol(xout$data)-2
  snames <- names(xout$data)
  xmad <- rep(NA,nsample)
  for(i in 2+(1:nsample)) {
    sout <- xout$output[xout$output$ID==snames[i],]
    xmad[i-2] <- mad(na.omit(xout$data[,i]) - rep(sout$seg.mean,sout$num.mark))
  }
  pfreq <- gain <- loss <- rep(0, nrow(xout$data))
  for(i in 1:nsample) {
#    ii <- !is.na(xout$data[,i+2])
    genomdat <- xout$data[,i+2]
# ii = location of the missing values and infinity
    ii <- which(is.finite(genomdat))
# segment means as a vector
    segout <- xout$output[xout$output$ID==snames[i+2],]
    segmean <- rep(segout$seg.mean, segout$num.mark)
# gains and losses
    pfreq[ii] <- pfreq[ii] + 1
    gain[ii] <- gain[ii] + 1*((segmean - median(segmean))/xmad[i] > threshold)
    loss[ii] <- loss[ii] - 1*((segmean - median(segmean))/xmad[i] < -threshold)
  }
  out <- list()
  out$chrom <- xout$data$chrom
  out$maploc <- xout$data$maploc
  out$pfreq <- pfreq
  out$gain <- gain/pfreq
  out$loss <- loss/pfreq
  as.data.frame(out)
}
