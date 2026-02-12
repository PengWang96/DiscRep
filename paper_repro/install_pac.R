args <- commandArgs(trailingOnly = TRUE)
quiet <- any(args == "--quiet")
pkg_name <- "DiscRep"

if (!requireNamespace("pkgbuild", quietly = TRUE)) {
  stop("Package 'pkgbuild' is required. Install it with install.packages('pkgbuild').")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required. Install it with install.packages('devtools').")
}

devtools::document()
# devtools::test()
devtools::check()
pkgbuild::clean_dll()

# Clear loaded namespace first; corrupt lazy-load db can break unload in devtools::install().
if (pkg_name %in% loadedNamespaces()) {
  try(unloadNamespace(pkg_name), silent = TRUE)
}

# Remove any existing installation folders from all library paths.
for (lib in .libPaths()) {
  pkg_dir <- file.path(lib, pkg_name)
  if (dir.exists(pkg_dir)) {
    unlink(pkg_dir, recursive = TRUE, force = TRUE)
  }
}

devtools::install(quiet = quiet)

library(DiscRep)
P_mis(0)
