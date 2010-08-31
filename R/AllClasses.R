setClassUnion("chromid", members =  c("numeric", "ordered"))
setClassUnion("index", members =  c("numeric", "logical", "character")) # from Matrix
setClassUnion("listOrNull", c("NULL","list"))

setClass("CNA", representation(chrom="chromid", maploc="numeric", genomdat="ANY", data.type="character", bigmemorydesc="listOrNull"))

setValidity("CNA", function(object) {
  msg <- NULL
  ismat <- is(slot(object, "genomdat"), "matrix")
  isbmat <- is(slot(object, "genomdat"), "big.matrix")
  isfmat <- is(slot(object, "genomdat"), "ff_matrix")
  if (!any(ismat,isbmat,isfmat)) {
    msg <- "genomdat has bad matrix type"
  } else {
    nc <- length(object@chrom)
    nm <- length(object@maploc)
    ng <- nrow(object@genomdat)
    if (nc != nm || nc != ng) 
      msg <- paste("\n  Unequal sizes: chrom = ", nc, "; maploc = ", nm, "; genomdat = ", ng, sep="")
  }
  if (is.null(msg)) TRUE else msg
})

setClass("DNAcopy", contains="CNA", representation(call="call", output="data.frame", segRows="data.frame"))
