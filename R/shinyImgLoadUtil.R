.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath(
    prefix = "assets",
    directoryPath = system.file(
      "assets",
      package = "CRISPRclassify"
    )
  )
}

.onUnload <- function(libname, pkgname) {
  shiny::removeResourcePath("assets")
}
