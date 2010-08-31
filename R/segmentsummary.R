segments.summary <- function(x) {
  if (class(x) != "DNAcopy") stop("First arg must be the result of segment")
  summary(x)
}
