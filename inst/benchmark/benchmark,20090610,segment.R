######################################################################
# Type: Redundancy test
# Created by: Henrik Bengtsson <hb@stat.berkeley.edu>
# Created on: 2009-06-10
######################################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Startup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
scriptName <- "benchmark,20090610,segment"

library("DNAcopy")
library("R.utils")

# Record current random seed
sample(1) # Assert that a random seed exists
oldSeed <- .Random.seed
# Alway use the same random seed
set.seed(0xbeef)

# Tolerance (maybe decrease?)
tol <- .Machine$double.eps^0.5

pd <- packageDescription("DNAcopy")
pkgStr <- sprintf("%s v%s", pd$Package, pd$Version) 

figPath <- Arguments$getWritablePath("figures")

benchmarkName <- paste(c(scriptName, gsub(" ", "_", pkgStr)), collapse=",")

logFilename <- sprintf("%s.log", benchmarkName)
log <- Verbose(logFilename, threshold=-10, timestamp=TRUE)

log && header(log, "BENCHMARKING")
log && cat(log, "Script: ", scriptName)
log && print(log, sessionInfo())

benchmarkFilename <- sprintf("%s.Rbin", benchmarkName)
force <- FALSE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Main benchmarking loop
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sizes of data sets to be benchmarked
Js <- c(1e3, 1e4, 1e5, 2e5, 5e5, 1e6)


if (!force && isFile(benchmarkFilename)) {
  benchmarkData <- loadObject(benchmarkFilename)
} else {
  benchmarkData <- data.frame(J=NULL, seg=NULL, weightSeg=NULL)
}


for (jj in seq(along=Js)) {
  # Number of loci
  J <- as.integer(Js[jj])

  log && enter(log, sprintf("Case #%d (J=%d) of %d", jj, J, length(Js)))

  if (is.element(J, benchmarkData$J)) {
    log && cat(log, "Already done.")
    log && exit(log)
    next
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Simulating copy-number data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  x <- sort(runif(J, min=0, max=1000))
  w <- runif(J)
  mu <- double(J)
  jj <- (200 <= x & x < 300)
  mu[jj] <- mu[jj] + 1
  jj <- (650 <= x & x < 800)
  mu[jj] <- mu[jj] - 1
  w[jj] <- 0.001
  eps <- rnorm(J, sd=1/2)
  y <- mu + eps

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up a raw CNA object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cnR <- CNA(
    genomdat  = y,
    chrom     = rep(1, times=J),
    maploc    = x,
    data.type = "logratio",
    sampleid  = "SampleA"
  )
  log && print(log, cnR)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Non-weighted segmentation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && enter(log, "Non-weighted segmentation")
  t1 <- system.time({
  fitR <- segment(cnR, verbose=1)
  })[3]
  log && printf(log, "Processing time: %.3f secs\n", t1)
  log && print(log, fitR)
  log && exit(log)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weighted segmentation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && enter(log, "Weighted segmentation")
  t2 <- system.time({
  fitR <- segment(cnR, weights=w, verbose=1)
  })[3]
  log && printf(log, "Processing time: %.3f secs\n", t1)
  log && print(log, fitR)
  log && exit(log)
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Record benchmarking
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  timings <- data.frame(J=J, seg=t1, weightSeg=t2)
  benchmarkData <- rbind(benchmarkData, timings)
  log && print(log, benchmarkData)

  # Saving to file  
  saveObject(benchmarkData, file=benchmarkFilename)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cleanup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reset to previous random seed
  .Random.seed <- oldSeed

  log && exit(log)
} # for (jj ...)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Benchmarking summary
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
log && print(log, benchmarkData)

log && header(log, "APPENDIX")
log && print(log, sessionInfo())


figName <- paste(c(scriptName, gsub(" ", "_", pkgStr)), collapse=",")
width <- 640
height <- 0.618*width
filename <- sprintf("%s.png", figName)
pathname <- file.path(figPath, filename)
devNew(png, pathname, width=width, height=height)
n <- ncol(benchmarkData)-1
matplot(benchmarkData[1], benchmarkData[,-1], type="b", pch=20, lwd=3, 
                               xlab="J", ylab="seconds", main=pkgStr)
legend("topleft", colnames(benchmarkData)[-1], col=1:n, lty=1:n, lwd=3)
devDone()


######################################################################
# HISTORY:
# 2009-06-10
# o Benchmarking show a major improvement in the algorithm when going
#   from DNAcopy v1.19.0 to the recent DNAcopy v1.19.2.  It was 
#   roughly O(J*ln(J)) and now it is O(J).  For a chromosome with
#   500,000 loci, we observed a speed up in the weighted case going
#   from 20 mins to 30 seconds, which is a 40 times speedup.
# o Created.
######################################################################
