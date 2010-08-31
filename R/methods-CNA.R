# accessor and replacement methods for CNA
setMethod("$", "CNA", function(x, name) {
  if (name == "chrom" | name == "maploc") {
    slot(x, name)
  } else {
    i4namedsample <- match(name, colnames(x@genomdat))
    if (is.na(i4namedsample)) {
      stop(paste("\"", name, "\" doesn't exist", sep=""))
    } else {
      x@genomdat[,i4namedsample]
    }
  }
})

setReplaceMethod("$", "CNA", function(x, name, value) { 
  if (name == "chrom" | name == "maploc") {
    slot(x, name) <- value
    x
  } else {
    stop("only chrom and maploc can be modified using the accessor $<-.")
  }
})

setMethod("[", signature(x="CNA", i="missing", j="missing", drop="ANY"),
          function(x, i, j, ..., drop) x@genomdat[])

setMethod("[", signature(x="CNA", i="index", j="missing", drop="missing"),
          function(x, i, j, ..., drop) x@genomdat[i,])

setMethod("[", signature(x="CNA", i="missing", j="index", drop="missing"),
          function(x, i, j, ..., drop) x@genomdat[,j])

setMethod("[", signature(x="CNA", i="index", j="index", drop="missing"),
          function(x, i, j, ..., drop) x@genomdat[i,j])

setMethod("[", signature(x="CNA", i="ANY", j="ANY", drop="ANY"),
          function(x, i, j, ..., drop) stop("drop is not valid for CNA"))

setReplaceMethod("[", signature(x="CNA", i="index", j="index"),
                 function(x, i, j, value) {x@genomdat[i,j] <- value; x})

setReplaceMethod("[", signature(x="CNA", i="missing", j="index"),
                 function(x, i, j, value) {x@genomdat[,j] <- value; x})

setReplaceMethod("[", signature(x="CNA", i="index", j="missing"),
                 function(x, i, j, value) {x@genomdat[i,] <- value; x})

setReplaceMethod("[", signature(x="CNA", i="missing", j="missing"),
                 function(x, i, j, value) {x@genomdat[] <- value; x})

# print and show methods
setMethod("print", "CNA", function(x, ...) {
  cat("CNA object of type:", attr(x,"data.type"),
      "\nNumber of Probes:  ", nrow(x@genomdat),
      "\nNumber of Samples: ", ncol(x@genomdat),"\n")
})

setMethod("show", "CNA", function(object){
  cat(attr(object,"data.type"), "CNA data with", ncol(object@genomdat), "samples and", nrow(object@genomdat), "probes \n") 
})

# function to initialize a CNA object
CNA <- function(genomdat, chrom, maploc, data.type=c("logratio","binary"),
                sampleid=NULL) {
  data.type <- match.arg(data.type)
  if (is(genomdat,"ff_matrix")) {
    bigmem <- FALSE
    fbdata <- TRUE
    bfile <- physical(genomdat)$filename
    bfile <- basename(bfile)
    if (!file.exists(bfile))
      stop(paste("data file not in current working directory:", getwd()))
  } else if (is(genomdat, "big.matrix")) {
    bigmem <- TRUE
    fbdata <- TRUE
    if (is.nil(genomdat@address))
      stop("pointer address missing for the backing file of genomdat")
    fbspecs <- describe(genomdat)
    bigmemdesc <- fbspecs@description
    bfile <- bigmemdesc$filename
    if (!file.exists(bfile))
      stop(paste("data file not in current working directory:", getwd()))
  } else {
    bigmem <- FALSE
    fbdata <- FALSE
    if (is.data.frame(genomdat)) genomdat <- as.matrix(genomdat)
    if (!is.numeric(genomdat)) stop("genomdat must be numeric")
    if (is.vector(genomdat)) genomdat <- as.matrix(genomdat)
  }
  nsample <- ncol(genomdat)
  if (is.factor(chrom) | is.character(chrom)) chrom <- orderedChrom(chrom)
  if (!is.numeric(maploc)) stop("maploc must be numeric")

  ina <- (!is.na(chrom) & is.finite(maploc))
  if (sum(!ina)>0)
    warning("markers with missing chrom and/or maploc removed\n")

  if (fbdata) {
    sortindex <- c(which(ina)[order(chrom[ina], maploc[ina])], which(!ina))
  } else {
    sortindex <- which(ina)[order(chrom[ina], maploc[ina])]
  }
  chrom <- chrom[sortindex]
  maploc <- maploc[sortindex]
  for (i in 1:nsample) genomdat[,i] <- genomdat[sortindex,i]
  
  if (!missing(sampleid)) {
    if (length(sampleid) != ncol(genomdat)) {
      warning("length(sampleid) and ncol(genomdat) differ, names ignored\n")
      sampleid <- paste("Sample", 1:ncol(genomdat))
    } 
  } else {
    cnames <- colnames(genomdat)
    if (!is.null(cnames)) {
      if (any(duplicated(cnames)) | any(cnames=="")) {
        warning("replacing ill defined colnames(genomdat)")
        sampleid <- paste("Sample", 1:ncol(genomdat))
      } else {
        message("using colnames(genomdat) as sampleid")
        sampleid <- cnames
      }
    } else {
      sampleid <- paste("Sample", 1:ncol(genomdat))
    }
  }

  if (is(genomdat, "big.matrix")) {
    bigmemdesc$colNames <- sampleid
  } else {
    colnames(genomdat) <- sampleid
  }

  if (length(ii <- which(diff(maploc)==0)) > 0) {
    if (any(chrom[ii]==chrom[ii+1])) warning("array has repeated maploc positions\n")
  }

  if (bigmem) {
    zzz <- new("CNA", chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type, bigmemorydesc=bigmemdesc)
  } else {
    zzz <- new("CNA", chrom=chrom, maploc=maploc, genomdat=genomdat, data.type=data.type)
  }
  zzz
}

# utility function to make non-numeric chrom as ordered data 
orderedChrom <- function(x) {
  if (is.factor(x)) ux <- levels(x)
  else ux <- unique(x)
  nux <- suppressWarnings(as.numeric(ux))
  ordered(x, levels=c(sort(nux[!is.na(nux)]),sort(ux[is.na(nux)])))
}

# subset method for CNA object
setMethod("subset", "CNA", function(x, chromlist=NULL, samplelist=NULL, fbdata=NULL, backingfile=NULL, ...) {
  if (!is(slot(x,"genomdat"), "matrix")) {
    if (missing(fbdata)) stop("genomdat is filebacked; use fbdata option")
    if (!is.logical(fbdata)) stop("fbdata should be of type logical")
    if (fbdata) {
      if (missing(backingfile) | is.null(backingfile)) stop("you must specify a backing file")
      if (!is.character(backingfile)) stop("backingfile should be a character string")
    }
  } else {
    if (missing(fbdata) | is.null(fbdata)) fbdata <- FALSE
  }
  
  uchrom <- unique(x@chrom)
  if (missing(chromlist) | is.null(chromlist)) chromlist <- uchrom
  if (length(setdiff(chromlist, uchrom)) > 0)
    stop("chromlist contains chromosomes not in the data")
  if (length(chromlist) > length(unique(chromlist)))
    warning("duplicate chromosomes in chromlist ignored")
  sampleid <- colnames(x@genomdat)
  if (missing(samplelist) | is.null(samplelist)) samplelist <- sampleid
  nsample <- length(sampleid)
  if (length(setdiff(samplelist, 1:nsample)) > 0 & length(setdiff(samplelist, sampleid)) > 0)
    stop("samplelist should be a list of valid sample numbers or names")
  if (!is.numeric(samplelist)) samplelist <- match(samplelist, sampleid)
  if (length(samplelist) > length(unique(samplelist)))
    warning("duplicate samples in samplelist removed")
  samplelist <- unique(samplelist)
  ii <- which(x@chrom %in% chromlist)

  if (fbdata) {
    nr <- length(ii)
    nc <- length(samplelist)
    if (is(slot(x,"genomdat"), "ff_matrix")) {
      genomdat <- ff(vmode="double", dim=c(nr, nc), file=backingfile)
      colnames(genomdat) <- sampleid[samplelist]
      for (i in 1:nc) genomdat[,i] <- x@genomdat[ii,samplelist[i]]
      zzz <- new("CNA", chrom=x@chrom[ii], maploc=x@maploc[ii], genomdat=genomdat, data.type=x@data.type)
    } else {
      descfile <- paste(backingfile, "bmdesc", sep=".")
      genomdat <- filebacked.big.matrix(nr, nc, backingfile=backingfile, descriptorfile=descfile, dimnames=list(c(), sampleid[samplelist]))
      bigmemdesc <- genomdat@description
      for (i in 1:nc) genomdat[,i] <- x@genomdat[ii,samplelist[i]]
      unlink(descfile)
      zzz <- new("CNA", chrom=x@chrom[ii], maploc=x@maploc[ii], genomdat=genomdat, data.type=x@data.type, bigmemorydesc=bigmemdesc)
    }
  } else {
    zzz <- new("CNA", chrom=x@chrom[ii], maploc=x@maploc[ii], genomdat=x@genomdat[ii,samplelist,drop=F], data.type=x@data.type)
  }
  zzz
})

# attach big.matrix if pointer is nil allows to save objects as Rdata
attach.genomdat <- function(obj) {
  if (is(obj@genomdat, "matrix")) {
    message("genomdat not file backed; nothing to attach")
  } else if (is(obj@genomdat, "big.matrix")) {
    if (is.nil(obj@genomdat@address))
      bmdesc <- new("big.matrix.descriptor", description=obj@bigmemorydesc)
      obj@genomdat <- attach.resource(bmdesc)
  } else {
    open.ff(obj@genomdat)
  }
}

# Chromosome.Lengths <- c(263, 255, 214, 203, 194, 183, 171, 155, 145, 144, 144, 143, 114, 109, 106, 98, 92, 85, 67, 72, 50, 56, 164, 59)
# names(Chromosome.Lengths) <- c(as.character(1:22),"X","Y")
