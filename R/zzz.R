.onLoad <- function(libname, pkgname) {
  library.dynam("DNAcopy", pkgname, libname)
  packageStartupMessage("\n**************************************************************************\n   The data format for the CNA object will be changed in version 1.23.0 \n     Instead of a data frame it will be a list of 3 (or more) objects \n     chrom and maploc will be vectors and CN/LOH data will be a matrix \n**************************************************************************\n")
}
