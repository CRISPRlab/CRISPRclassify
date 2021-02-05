.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath(
    prefix = "csvRef",
    directoryPath = system.file(
      "csvRef",
      package = "CRISPRclassify"
    )
  )

  shiny::addResourcePath(
    prefix = "assets",
    directoryPath = system.file(
      "assets",
      package = "CRISPRclassify"
    )
  )

  shiny::addResourcePath(
    prefix = "lib",
    directoryPath = system.file(
      "lib",
      package = "CRISPRclassify"
    )
  )

  shiny::addResourcePath(
    prefix = "models",
    directoryPath = system.file(
      "models",
      package = "CRISPRclassify"
    )
  )

}

.onUnload <- function(libname, pkgname) {
  shiny::removeResourcePath("csvRef")
  shiny::removeResourcePath("assets")
  shiny::removeResourcePath("lib")
  shiny::removeResourcePath("models")
}
