rm(list = ls())
library(DiscRep)
library(ggplot2)
library(foreach)
library(doParallel)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

t1 <- Sys.time()
set.seed(234)

# ---- Common parameters ----
k <- 0.25
bbar <- 0.1
num_simulations <- 2000
rnse <- 1 # replication noise level/se
sample_size <- 100
m <- 10
bc_grp_num <- 4 # batch contaminated groups

# ---- k_vec grid (shared across all bb_sd) ----
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate

# ---- Loop over bb_sd values ----
bb_sd_vec <- c(0, 0.4, 0.8, 1.2)

for (bb_sd in bb_sd_vec) {
  cat("\n========== bb_sd =", bb_sd, "==========\n")

  # --- Simulate batch data ---
  Gv <- t(sapply(1:m, function(x) rbinom(sample_size, 1, 0.4)))
  Batch <- t(apply(Gv, 1, function(x) shuffle(x, rep = floor(sample_size * 0.2))))
  rst <- t(sapply(
    1:num_simulations,
    function(x) {
      sim_batch(
        bb_sd = bb_sd, k = k, bbar = bbar, m,
        bc_grp_num, sample_size, Gv, Batch, rnse
      )
    }
  ))
  nrst <- ncol(rst)

  hat_beta <- rst[, seq(1, nrst, 2)]
  hat_sigma_sq <- rst[, seq(2, nrst, 2)]^2

  saveRDS(list(hat_beta = hat_beta, hat_sigma_sq = hat_sigma_sq),
    file = paste0("../output/hat_beta_sigma_sq_eta_", bb_sd, ".rds")
  )

  # --- MCMC (parallel) ---
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  sim_results <- foreach(x = 1:num_simulations, .combine = combine_results, .packages = "DiscRep") %dopar% {
    metropolis_hastings(N, r, m, hat_beta[x, ], hat_sigma_sq[x, ], k_vec = k_vec)
  }
  stopCluster(cl)

  # --- Extract results ---
  p_values <- sim_results$p_value
  bar_betas <- sim_results$bar_beta
  sen <- sum(p_values <= 0.05) / length(p_values)

  saveRDS(sim_results,
    file = paste0("../output/sim_results_eta_", bb_sd, ".rds")
  )

  # --- Summary ---
  cat("  Sensitivity:", sen, "\n")
  cat("  Acceptance rate:", sim_results[["acception"]], "\n")
  cat("  p-value range: [", min(p_values), ",", max(p_values), "]\n")
  cat("  p-value mean:", mean(p_values), "\n")
}

t2 <- Sys.time()
cat("\nTotal elapsed time:", format(t2 - t1), "\n")
