plotSample <- function(x, sampleid=NULL, chromlist=NULL, xmaploc=FALSE,
                       col=c("black","green"), pch=".", cex=NULL, altcol=TRUE,
                       segcol="red", lwd=3, zeroline=TRUE, zlcol="grey",
                       xlab=NULL, ylab=NULL, main=NULL, ...) {
  if (class(x) != "DNAcopy") stop("First arg must be a DNAcopy object")
  if (missing(sampleid)) {sampleid <- 1}
  subx <- subset(x, chromlist=chromlist, samplelist=sampleid[1])
# get the data for plotting
  genomdat <- subx$data[,3]
  ina <- is.finite(genomdat)
  genomdat <- genomdat[ina]
  chrom <- subx$data[ina,1]
  uchrom <- unique(chrom)
  segres <- subx$output
# setup the X-axis based on xmaploc
  if (xmaploc) {
    maploc <- subx$data[ina,2]
    rmaploc <- sapply(uchrom, function(i, maploc, chrom) range(maploc[chrom==i]), maploc, chrom)
    nc <- length(uchrom)
    if ((nc>1) && any(rmaploc[1,-1] < rmaploc[2,-nc])) {
      cmaploc <- cumsum(as.numeric(rmaploc[2,]))
      for (i in 2:nc) {
        maploc[chrom==uchrom[i]] <- cmaploc[i-1] + maploc[chrom==uchrom[i]] 
      }
    }
    xlabel <- "Genomic Position"
  } else {
    maploc <- 1:sum(ina)
    xlabel <- "Index"
  }
# setup altenating colors
  if (altcol & length(uchrom)>1) {
    colvec <- rep(1, length(chrom))
    j <- 0
    for(i in uchrom) {
      j <- (j+1) %% 2
      colvec[chrom == i] <- j+1
    }
  } else {
    colvec <- 1
  }
# set other graphical parameters
  if (missing(cex)) cex <- ifelse(pch == ".", 3, 1)
  if (missing(main)) main <- names(subx$data)[3]
  if (missing(xlab)) xlab <- xlabel
  if (missing(ylab)) {
    if (attr(subx$data, "data.type") == "logratio") {ylab <- "log(relative CN)"}
    else {ylab <- "LOH"}
  }
# plot the data
  plot(maploc, genomdat, col=col[colvec], pch=pch, cex=cex, main=main, xlab=xlab, ylab=ylab, ...)
# add the segment means
  ii <- cumsum(c(0, segres$num.mark))
  mm <- segres$seg.mean
  kk <- length(ii)
  segments(maploc[ii[-kk]+1], segres$seg.mean, x1=maploc[ii[-1]], y1=segres$seg.mean, col = segcol, lwd=lwd)
#  for (i in 1:(kk - 1)) {
#    lines(maploc[c(ii[i]+1,ii[i+1])], rep(mm[i], 2), col = segcol, lwd=lwd)
#  }
# add the zeroline
  if (zeroline) abline(h=0, col=zlcol, lwd=lwd)
}
