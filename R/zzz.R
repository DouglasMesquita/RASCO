.onLoad <- function(libname, pkgname) {
  ##-- Installing INLA ifnecessary
  if(!requireNamespace("INLA", quietly = TRUE)) {
    utils::install.packages("INLA",
                     repos = c(getOption("repos"),
                               INLA = "https://inla.r-inla-download.org/R/stable"),
                     dep = TRUE)
  }

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
