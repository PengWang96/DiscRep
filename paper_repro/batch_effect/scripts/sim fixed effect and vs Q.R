# k = 0, fixed effect, compare posterior-PRP and Q test 
# based on batch effect simulation
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
library(metafor)
bb_sd_vec <- c(0, 0.4, 0.8)

k <- 0
bbar <- 0.1
num_simulations <- 2000 #
rnse <- 1.0 # replication noise level/se
sample_size <- 100
m <- 10
bc_grp_num <- 4 # batch contaminated groups
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate

CochranQ <- function(hat_beta, hat_sigma_sq) {
  res <- rma(yi = hat_beta, sei = sqrt(hat_sigma_sq), method = "FE")
  p_values <- res$QEp
  return(p_values)
}

numCores <- detectCores()
pvalue_results <- list()

for (bb_sd in bb_sd_vec) {
  Gv <- t(sapply(1:m, function(x) rbinom(sample_size, 1, 0.4)))
  Batch <- t(apply(Gv, 1, function(x) shuffle(x, rep = floor(sample_size * 0.20))))
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

  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  sim_results <- foreach(x = 1:num_simulations, .combine = combine_results, .packages = "DiscRep") %dopar% {
    metropolis_hastings(N, r, m, hat_beta[x, ], hat_sigma_sq[x, ], k_vec = k_vec)
  }
  stopCluster(cl)
  p_values_prp <- sim_results$p_value

  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  p_values_q <- foreach(x = 1:num_simulations, .combine = c, .packages = "metafor") %dopar% {
    CochranQ(hat_beta[x, ], hat_sigma_sq[x, ])
  }
  stopCluster(cl)

  bb_sd_tag <- gsub("\\.", "_", as.character(bb_sd))
  pvalue_results[[paste0("PRP_bb_sd_", bb_sd_tag)]] <- p_values_prp
  pvalue_results[[paste0("Q_bb_sd_", bb_sd_tag)]] <- p_values_q
}

saveRDS(pvalue_results, file = "../output/pvalues_k_0.rds")
t2 <- Sys.time()
(t2 - t1)