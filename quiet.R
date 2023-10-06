# quiet function
### dependency for ewas_qgcomp function, quiets summary output
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
