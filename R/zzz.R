.onLoad <- function(libname, pkgname) {
  INLA:::inla.dynload.workaround()
}
