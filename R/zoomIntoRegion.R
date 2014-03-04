zoomIntoRegion <- function(x, chrom, sampleid, maploc.start=NULL, maploc.end=NULL, pt.pch=NULL, pt.cex=NULL, pt.col=NULL, segcol=NULL, seglwd=NULL, main=NULL, xlab=NULL, ylab=NULL, ...) {
  if (class(x) != "DNAcopy") stop("First arg must be a DNAcopy object")
  tmp <- subset(x, chrom=chrom[1], samplelist=sampleid[1])
  lrdata <- tmp$data
  if (missing(maploc.start)) maploc.start <- min(lrdata$maploc, na.rm=T) - 1
  if (missing(maploc.end)) maploc.end <- max(lrdata$maploc, na.rm=T) + 1
  ii <- ((lrdata$maploc >= maploc.start) & (lrdata$maploc <= maploc.end))
  if (missing(pt.pch)) pt.pch <- "."
  if (missing(pt.cex))
      pt.cex <- ifelse(pt.pch==".", 3, 1)
  if (missing(pt.col)) pt.col <- "green3"
  if (missing(segcol)) segcol <- "red"
  if (missing(seglwd)) seglwd <- 3
  if (missing(main)) 
      main <- paste("chr", chrom, ": ", maploc.start,"-", maploc.end, " from sample ", sampleid, sep="")
  if (missing(xlab)) xlab = "Genomic Position"
  if (missing(ylab)) ylab = "log-ratio"
  
  plot(lrdata[ii,2], lrdata[ii,3], main = main, xlab=xlab, ylab = ylab, pch = pt.pch, cex = pt.cex, col = pt.col, ...)
  segs <- tmp$output
  jj <- ((segs$loc.start <= maploc.end) & (segs$loc.end >= maploc.start))
  segs <- segs[jj,]
  k <- nrow(segs)
  segs$loc.start[1] <- maploc.start
  segs$loc.end[k] <- maploc.end
  segments(segs$loc.start, segs$seg.mean, x1=segs$loc.end, y1=segs$seg.mean, col = segcol, lwd = seglwd)
#  for(i in 1:k) {
#    lines(c(segs$loc.start[i],segs$loc.end[i]), rep(segs$seg.mean[i],2), col=segcol, lwd=seglwd)
#  }
}
