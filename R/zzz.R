.onLoad <- function(libname, pkgname) {
  ##-- INLA workaround for xenial
  if(Sys.info()['sysname'] == "Linux") {
    release <- system(command = "lsb_release -a", intern = TRUE, ignore.stderr = TRUE)
    release <- paste(gsub(x = release, pattern = "\t", replacement = " "), collapse = " ")

    ubuntu_version <- gsub(x = release, pattern = ".*Codename: ([a-zA-Z]*)$", replacement = "\\1")

    if(ubuntu_version == "xenial") {
      INLA:::inla.dynload.workaround()
    }
  }
}
