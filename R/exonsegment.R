exon.segment <- function(gene, eloc, edat, ngrid=100, tol=1e-6) {
  ii <- order(gene, eloc)
  gene <- gene[ii]
  eloc <- eloc[ii]
  if (is.matrix(edat)) {
    edat <- edat[ii,]
  } else {
    edat <- cbind(edat[ii])
  }
  ugene <- unique(gene)
  ngene <- length(ugene)
  nsample <- ncol(edat)
  out.stat <- out.loc <- out.p <- matrix(0, ngene, nsample)
  ss <- 3*(1:nsample)
  for(i in 1:ngene) {
    exondat <- edat[gene==ugene[i],]
    gout <- exon.changepoint(exondat, ngrid, tol)
    out.stat[i,] <- gout[[1]]
    out.loc[i,] <- gout[[2]]
    out.p[i,] <- gout[[3]]
  }
  rownames(out.stat) <- rownames(out.loc) <- rownames(out.p) <- ugene
  list(statistic=out.stat, location=out.loc, p.value=out.p)
}

exon.changepoint <- function(exondat, ngrid=100, tol=1e-6) {
# 
# exondat -- is a matrix of normalized expression values
#   rows are ordered by location and columns are samples
#
  nsample <- ncol(exondat) # number of samples
  n <- nrow(exondat) # number of exons in the gene
# initialize sample specific output
  estat <- epval <- eloc <- rep(0, nsample)
# calculate the max t-stat, location and p-value
  for(i in 1:nsample) {
    exondati <- exondat[,i]
#   center the data
    exondati <- (exondati - mean(exondati))
# call the p-value subroutine
    zzz <- .Fortran("esegp",
                    as.integer(n),
                    as.double(exondati),
                    ostat=double(1),
                    eloc=integer(1),
                    pval=double(1),
                    as.integer(ngrid),
                    as.double(tol),
                    PACKAGE="DNAcopy")
    estat[i] <- zzz$ostat
    epval[i] <- zzz$pval
    eloc[i] <- zzz$eloc
  }
  list(estat, eloc, epval)
}
