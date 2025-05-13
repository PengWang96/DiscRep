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
bb_sd <- 0.8 # 0 # 0.4 # 

k <- 0
bbar <- 0.1
num_simulations <- 2000 #
rnse <- 1.0  # replication noise level/se
sample_size <- 100
m <- 10
bc_grp_num = 4 # batch contaminated groups
Gv = t(sapply(1:m, function(x) rbinom(sample_size, 1, 0.4)))
Batch = t(apply(Gv,1, function(x) shuffle(x, rep=floor(sample_size*0.20))))
rst <- t(sapply(1:num_simulations, 
                function(x) sim_batch(bb_sd = bb_sd, k = k, bbar = bbar, m, 
                                      bc_grp_num, sample_size, Gv, Batch, rnse)))
nrst <- ncol(rst)

hat_beta <- rst[, seq(1, nrst, 2)]
hat_sigma_sq <- rst[, seq(2, nrst, 2)]^2

pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)
sim_results <- foreach(x = 1:num_simulations, .combine = combine_results, .packages = "DiscRep") %dopar% {
  metropolis_hastings(N, r, m, hat_beta[x,], hat_sigma_sq[x,], k_vec = k_vec)
}
stopCluster(cl)


# get p-value and bar_beta
p_values <- sim_results$p_value
bar_betas <- sim_results$bar_beta


save(p_values,
     file = paste0("../output/pvalue_eta_", bb_sd, "_k_0.rda"))
ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(-0.00, 1.00)) +
  theme_bw() +
  theme(plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion") +
  ylab("Count")

sim_results[["acception"]]
mean(p_values)
sum(p_values <= 0.05)/length(p_values)
t2 <- Sys.time()
(t2-t1)







rm(list = setdiff(ls(), c("hat_beta", "hat_sigma_sq", "bb_sd")))
t1 <- Sys.time()
library(metafor)
set.seed(234)

m <- ncol(hat_beta)
num_simulations <- nrow(hat_beta)

boxplot(hat_beta)

CochranQ <- function(hat_beta, hat_sigma_sq) {
  res <- rma(yi = hat_beta, sei = sqrt(hat_sigma_sq), method = "FE")
  p_values <- res$QEp
  return(p_values)
}
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)
p_values <- foreach(x = 1:num_simulations, .combine = c, .packages = "metafor") %dopar% {
  CochranQ(hat_beta[x,], hat_sigma_sq[x,])
}
stopCluster(cl)
save(p_values,
     file = paste0("../output/Frequentist_pvalue_eta_", bb_sd, ".rda"))


ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(0.00, 1.00)) +
  theme_bw() +
  theme(plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Frequentist p-values")

t2 <- Sys.time()
(t2-t1)




# bb_sd <- 0 # 0.4 # 0.8
rm(list = setdiff(ls(), "bb_sd"))
library(ggplot2)
load(paste0("../output/Frequentist_pvalue_eta_", bb_sd, ".rda"))
Q_p_values <- p_values
sum(Q_p_values <= 0.1)/length(Q_p_values)
load(paste0("../output/pvalue_eta_", bb_sd, "_k_0.rda"))
sum(p_values <= 0.1)/length(p_values)
dataa <- data.frame(Q = Q_p_values, PRP = p_values)
plot <- ggplot(dataa, aes(x = Q, y = PRP)) +
  geom_point(color = "#6186ad") +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, color = "black") +
  xlab("Cochran's Q P-values") + 
  ylab("Posterior-PRPs") + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 20), # Increase x-axis title font size
    axis.title.y = element_text(size = 20)  # Increase y-axis title font size
  )
plot
ggsave(paste0("../plot/fixedPRPvsQ_", bb_sd, ".pdf"), plot)
sum(p_values <= 0.1)/length(p_values)
sum(Q_p_values <= 0.1)/length(Q_p_values)

