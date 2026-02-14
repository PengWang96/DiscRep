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
bb_sd <- 0.4 # 1.2 # 0 # 0.8 # 

k <- 0.25
bbar <- 0.1
num_simulations <- 2000
rnse <- 1  # replication noise level/se
sample_size <- 100
m <- 10
bc_grp_num <- 4 # batch contaminated groups
Gv <- t(sapply(1:m, function(x) rbinom(sample_size, 1, 0.4)))
Batch <- t(apply(Gv,1, function(x) shuffle(x, rep=floor(sample_size*0.2))))
rst <- t(sapply(1:num_simulations, 
                function(x) sim_batch(bb_sd = bb_sd, k = k, bbar = bbar, m, 
                                      bc_grp_num, sample_size, Gv, Batch, rnse)))
nrst <- ncol(rst)

hat_beta <- rst[, seq(1, nrst, 2)]
hat_sigma_sq <- rst[, seq(2, nrst, 2)]^2
boxplot(hat_beta)


save(hat_beta,
     file = paste0("../output/hat_beta_eta_", bb_sd, ".rda"))
save(hat_sigma_sq,
     file = paste0("../output/hat_sigma_sq_eta_", bb_sd, ".rda"))


# pvec <- seq(0, 0.05, length = 1000)
# k_vec <- sapply(pvec, inverse_P_mis)
pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)
sim_results <- foreach(x = 1:num_simulations, .combine = combine_results, .packages = "DiscRep") %dopar% {
  metropolis_hastings(N, r, m, hat_beta[x,], hat_sigma_sq[x,], k_vec = k_vec)
  # metropolis_hastings(N, r, m, hat_beta[x,], hat_sigma_sq[x,], k_vec = 0.27)
} 
stopCluster(cl)


# get p-value and bar_beta
p_values <- sim_results$p_value
bar_betas <- sim_results$bar_beta
sen <- sum(p_values <= 0.05)/length(p_values)


save(p_values,
     file = paste0("../output/pvalue_eta_", bb_sd, ".rda"))


ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=40, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(-0.00, 1.00)) +
  theme_bw() +
  theme(plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion") +
  ylab("Count")

sen
sim_results[["acception"]]
max(p_values)
min(p_values)
mean(p_values)
t2 <- Sys.time()
(t2-t1)
