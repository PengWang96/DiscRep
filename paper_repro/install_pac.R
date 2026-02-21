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
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)
hist(k_vec)
hist(pvec)

pvec <- seq(0,  0.05, 0.0001)
k_vec <- sapply(pvec, inverse_P_mis)
hist(k_vec)
hist(pvec)

data("dat.slf")
m <- nrow(dat.slf)
hat_beta <- dat.slf$y
hat_sigma_sq <- dat.slf$s2

results_random_log_uniform_k <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq)

results_random_fixed_k <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, k_vec = 0.2726814)

results_random_uniform_pmis <- metropolis_hastings(
  10000, 0.05, m, hat_beta, hat_sigma_sq,
  k_vec_dist = "uniform_pmis"
)

results_random_trunc_beta <- metropolis_hastings(
  10000, 0.05, m, hat_beta, hat_sigma_sq,
  k_vec_dist = "trunc_beta_pmis",
  beta_shape1 = 2,
  beta_shape2 = 8
)

results_random_constant_pmis <- metropolis_hastings(
  10000, 0.05, m, hat_beta, hat_sigma_sq,
  k_vec_dist = "constant_pmis",
  p_mis_constant = 0.02
)

results_random_uniform_k <- metropolis_hastings(
  10000, 0.05, m, hat_beta, hat_sigma_sq,
  k_vec_dist = "uniform_k"
)
