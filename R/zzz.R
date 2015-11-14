#formatR::tidy_dir("R", indent = getOption("formatR.indent", 2), arrow = getOption("formatR.arrow", T))
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to mixpack. mixpack is a package for mixture component visualization")
  
}