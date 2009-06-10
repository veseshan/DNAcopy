######################################################################
# Type: Redundancy test
# Created by: Henrik Bengtsson <hb@stat.berkeley.edu>
# Created on: 2009-06-10
######################################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Startup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("DNAcopy")

# Record current random seed
sample(1) # Assert that a random seed exists
oldSeed <- .Random.seed
# Alway use the same random seed
set.seed(0xbeef)

# Tolerance (maybe decrease?)
tol <- .Machine$double.eps^0.5


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y)))
w <- runif(J)
w[650:800] <- 0.001


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
print(cnR)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test: Non-weighted segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitR <- segment(cnR, verbose=1)
print(fitR)

# Expected results [from dput(fitR$output)]
truth <- structure(list(ID=c("SampleA", "SampleA", "SampleA", "SampleA",
"SampleA"), chrom=c(1, 1, 1, 1, 1), loc.start=c(0.551678240299225,
207.80026470311, 293.691275408491, 659.097356721759, 814.257808262482
), loc.end=c(207.684752298519, 292.710821842775, 658.396144397557,
812.704775715247, 999.108265852556), num.mark=c(201, 99, 349,
151, 200), seg.mean=c(0.0164, 1.0474, -0.0227, -1.0813, -0.0612
)), .Names=c("ID", "chrom", "loc.start", "loc.end", "num.mark",
"seg.mean"), row.names=c(NA, -5L), class="data.frame")

stopifnot(all.equal(fitR$output, truth, tolerance=tol))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test: Weighted segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitR <- segment(cnR, weights=w, verbose=1)
print(fitR)

# Expected results [from dput(fitR$output)]
truth <- structure(list(ID=c("SampleA", "SampleA", "SampleA"), chrom=c(1,
1, 1), loc.start=c(0.551678240299225, 206.219613086432, 293.691275408491
), loc.end=c(205.489653628320, 292.710821842775, 999.108265852556
), num.mark=c(198, 102, 700), seg.mean=c(-0.0045, 1.0659,
-0.0425)), .Names=c("ID", "chrom", "loc.start", "loc.end",
"num.mark", "seg.mean"), row.names=c(NA, -3L), class="data.frame")

stopifnot(all.equal(fitR$output, truth, tolerance=tol))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Cleanup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reset to previous random seed
.Random.seed <- oldSeed
