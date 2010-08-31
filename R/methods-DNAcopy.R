# accessor and replacement methods for CNA
setMethod("$", "DNAcopy", function(x, name) {
  if (name == "output" | name == "segRows") {
    slot(x, name)
  } else {
    callNextMethod()
  }
})

setMethod("print", "DNAcopy", function(x, ...) {
  show(x)
  print(x@call)
  showSegRows <- list(...)[['showSegRows']]
  if (is.null(showSegRows)) showSegRows <- FALSE
  if (showSegRows) output <- cbind(x@output, x@segRows)
  else output <- x@output
  fullOutput <- list(...)[['fullOutput']]
  if (is.null(fullOutput)) fullOutput <- FALSE
  if (fullOutput) print(output)
  else head(output)
})

setMethod("show", "DNAcopy", function(object){
  cat("segmented", attr(object,"data.type"), "CNA data with",
      ncol(object@genomdat), "samples and", nrow(object@genomdat), "probes \n")
})

# subset method for DNAcopy
setMethod("subset", "DNAcopy", function(x, chromlist=NULL, samplelist=NULL, fbdata=NULL, backingfile=NULL, ...) {
  zzz0 <- callNextMethod(x, chromlist, samplelist, fbdata, backingfile)
  call <- match.call()
  call[[1]] <- as.name("subset")
  zzz <- new("DNAcopy", call=call, zzz0)

  sampleid <- colnames(x@genomdat)
  if (missing(samplelist)) samplelist <- sampleid
  samplelist <- unique(samplelist)
  if (is.numeric(samplelist)) samplelist <- sampleid[samplelist]

  oID <- x@output$ID
  oChrom <- x@output$chrom

  # first get the rows in the output corresponding to the samples
  # use c to coerce into vector in case sapply returns matrix instead of list
  srows <- c(unlist(sapply(samplelist, function(i, id) {which(id==i)}, oID)))
  # now subset the chromosomes
  uchrom <- unique(x@chrom)
  if (missing(chromlist) | is.null(chromlist)) chromlist <- uchrom
  scrows <- srows[oChrom[srows] %in% chromlist]

  zzz@output <- x@output[scrows,]
  rownames(zzz@output) <- NULL
  zzz@segRows <- x@segRows[scrows,]
  rownames(zzz@segRows) <- NULL
  zzz
})

# summary method to report median, mad and sd
setMethod("summary", "DNAcopy", function(object, ...) {
    xout <- object$output
    nsample <- ncol(object@genomdat)
    sampleid <- colnames(object@genomdat)
    seg.median <- seg.sd <- seg.mad <- rep(NA, nrow(xout))
    ll <- 0
    for (isamp in 1:nsample) {
# genomdat = logratio data of sample isamp
      genomdat <- object[, isamp]
# ina = location of the missing values and infinity
      ina <- which(is.finite(genomdat))
# subset out the missing & infinity locations
      genomdat <- genomdat[ina]
      seglen <- xout$num.mark[xout$ID == sampleid[isamp]]
      kk <- length(seglen)
      seg.sd[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), sd)
      seg.median[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), median)
      seg.mad[ll+(1:kk)] <- tapply(genomdat, rep(1:kk,seglen), mad)
      ll <- ll + kk
    }
    xout$seg.sd <- round(seg.sd, 4)
    xout$seg.median <- round(seg.median, 4)
    xout$seg.mad <- round(seg.mad, 4)
    xout
})

# function to plot a DNAcopy object
setMethod("plot", signature(x="DNAcopy"),
  function(x, sampleid=NULL, chromlist=NULL, xmaploc=FALSE,
           col=c("black","green"), pch=".", cex=1, altcol=TRUE,
           segcol="red", lwd=3, zeroline=TRUE, zlcol="grey",
           xlab="", ylab="", main="", ...) {
  if (missing(sampleid)) {sampleid <- 1}
  subx <- subset(x, chromlist=chromlist, samplelist=sampleid[1], fbdata=FALSE)
# get the data for plotting
  genomdat <- subx@genomdat[,1]
  ina <- is.finite(genomdat)
  genomdat <- genomdat[ina]
  chrom <- subx$chrom[ina]
  uchrom <- unique(chrom)
  segres <- subx$output
# setup the X-axis based on xmaploc
  if (xmaploc) {
    maploc <- subx$maploc[ina]
    rmaploc <- sapply(uchrom, function(i, maploc, chrom) range(maploc[chrom==i]), maploc, chrom)
    nc <- length(uchrom)
    if ((nc>1) && any(rmaploc[1,-1] < rmaploc[2,-nc])) {
      cmaploc <- cumsum(rmaploc[2,])
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
  if (pch == ".") cex <- 3
  if (main=="") main <- colnames(subx@genomdat)
  if (xlab=="") xlab <- xlabel
  if (ylab=="") {
    if (x@data.type=="logratio") {ylab <- "log(CN)"}
    if (x@data.type=="binary") {ylab <- "LOH"}
  }
# plot the data
  plot(maploc, genomdat[ina], col=col[colvec], pch=pch, cex=cex, main=main, xlab=xlab, ylab=ylab, ...)
# add the segment means
  ii <- cumsum(c(0, segres$num.mark))
  mm <- segres$seg.mean
  kk <- length(ii)
  for (i in 1:(kk - 1)) {
    lines(maploc[c(ii[i]+1,ii[i+1])], rep(mm[i], 2), col = segcol, lwd=lwd)
  }
# add the zeroline
  if (zeroline) abline(h=0, col=zlcol, lwd=lwd)
})
