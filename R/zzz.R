.onLoad <- function(libname, pkgname) {
  library.dynam("DNAcopy", pkgname, libname)
  packageStartupMessage("\n**************************************************************************\n   The plan to change the data format for CNA object has been postponed   \n in order to ensure backward compatibility with older versions of DNAcopy \n**************************************************************************\n")
}
