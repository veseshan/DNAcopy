zoomIntoRegion <- function(x, chrom, sampleid, maploc.start=NULL, maploc.end=NULL, pt.pch=NULL, pt.cex=NULL, pt.col=NULL, segcol=NULL, seglwd=NULL, ...) {
  if (class(x) != "DNAcopy") stop("First arg must be a DNAcopy object")
  subx <- subset(x, chromlist=chrom[1], samplelist=sampleid[1], fbdata=FALSE)
  maploc <- subx$maploc
  lrdata <- subx[,1]
  if (missing(maploc.start)) maploc.start <- min(maploc, na.rm=T) - 1
  if (missing(maploc.end)) maploc.end <- max(maploc, na.rm=T) + 1
  ii <- ((maploc >= maploc.start) & (maploc <= maploc.end))
  if (missing(pt.pch)) pt.pch <- "."
  if (missing(pt.cex)) {
    if (pt.pch==".") { pt.cex <- 3}
    else {pt.cex <- 1}
  }
  if (missing(pt.col)) pt.col <- "green3"
  if (missing(segcol)) segcol <- "red"
  if (missing(seglwd)) seglwd <- 3
  
  plot(maploc[ii], lrdata[ii], main = sampleid, xlab="Genomic Position", ylab = "log-ratio", pch = pt.pch, cex = pt.cex, col = pt.col, ...)
  segs <- subx$output
  jj <- ((segs$loc.start <= maploc.end) & (segs$loc.end >= maploc.start))
  segs <- segs[jj,]
  k <- nrow(segs)
  segs$loc.start[1] <- maploc.start
  segs$loc.end[k] <- maploc.end
  for(i in 1:k) {
    lines(c(segs$loc.start[i],segs$loc.end[i]), rep(segs$seg.mean[i],2), col=segcol, lwd=seglwd)
  }
}
