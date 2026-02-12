args <- commandArgs(trailingOnly = TRUE)
quiet <- any(args == "--quiet")

if (!requireNamespace("pkgbuild", quietly = TRUE)) {
  stop("Package 'pkgbuild' is required. Install it with install.packages('pkgbuild').")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required. Install it with install.packages('devtools').")
}

devtools::document()
devtools::test()
devtools::check()
pkgbuild::clean_dll()
remove.packages("DiscRep")
devtools::install()
devtools::install(quiet = quiet)

unlink(file.path(.libPaths()[1], "DiscRep"), recursive = TRUE, force = TRUE)
