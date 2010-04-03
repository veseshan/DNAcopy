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

print(sessionInfo())


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 1000

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
print(cnR)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test: Non-weighted segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
t <- system.time({
fitR <- segment(cnR, verbose=1)
})
cat("Processing time:\n")
print(t)
print(fitR)

# Expected results
# These were obtained by dput(fitR$output) using DNAcopy v1.19.0
truth <- structure(list(ID = c("SampleA", "SampleA", "SampleA", "SampleA",
"SampleA"), chrom = c(1, 1, 1, 1, 1), loc.start = c(1.36857712641358,
201.604291098192, 303.775111911818, 650.741211604327, 800.302447052673
), loc.end = c(199.083976913244, 301.066882908344, 647.42697100155,
798.971758922562, 999.329038895667), num.mark = c(209, 105, 337,
138, 211), seg.mean = c(0.0256, 1.0099, -0.0084, -0.9792, -0.0289
)), .Names = c("ID", "chrom", "loc.start", "loc.end", "num.mark",
"seg.mean"), row.names = c(NA, -5L), class = "data.frame")

stopifnot(all.equal(fitR$output, truth, tolerance=tol))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test: Weighted segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
t <- system.time({
fitR <- segment(cnR, weights=w, verbose=1)
})
cat("Processing time:\n")
print(t)
print(fitR)

# Expected results
# These were obtained by dput(fitR$output) using DNAcopy v1.19.0
truth <- structure(list(ID = c("SampleA", "SampleA", "SampleA"), chrom = c(1,
1, 1), loc.start = c(1.36857712641358, 201.604291098192, 303.775111911818
), loc.end = c(199.083976913244, 301.066882908344, 999.329038895667
), num.mark = c(209, 105, 686), seg.mean = c(0.0259, 1.0004,
-0.0233)), .Names = c("ID", "chrom", "loc.start", "loc.end",
"num.mark", "seg.mean"), row.names = c(NA, -3L), class = "data.frame")

stopifnot(all.equal(fitR$output, truth, tolerance=tol))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Cleanup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reset to previous random seed
.Random.seed <- oldSeed

print(sessionInfo())


######################################################################
# HISTORY
# 2009-06-10
# o ROBUSTNESS: Added this test to assert that DNAcopy v1.19.2 and 
#   newer will numerically give the same results as DNAcopy v1.19.0.
#   This test is ran each time with R CMD check.
# o Created.
######################################################################
