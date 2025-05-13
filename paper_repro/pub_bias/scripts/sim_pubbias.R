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
set.seed(1234)
num_simulations <- 2000
c <- 20 # 0 # 10 # 20
k <- 0.25
# #### 20 studies
# #### sample size in each study
studynum <- c(10, 6, 4)
indnum <- c(50, 100, 120)
m <- sum(studynum)

samplesize <- unlist(sapply(1:length(studynum), function(x) rep(indnum[x], studynum[x])))

bbar <- log(2/3)
# bbar <- log(1/4)
beta <- bbar + rnorm(length(samplesize), sd = sqrt(k^2*bbar^2))
or <- exp(beta)

betafinal <- c()
sdfinal <- c()
pubfinal <- c()

pubbias <- matrix(0, num_simulations, m)
for (jj in 1:num_simulations){

  betalist<-c()
  sdlist<-c()
  for (j in 1:length(samplesize)){
    curor<-or[j]

    cur_indnum0 <- (samplesize[j] * 2) * 0.6
    cur_indnum1 <- (samplesize[j] * 2) * 0.4
    beta_cur_est<-c()
    sd_cur_est<-c()

    while (TRUE){
      pcontrol <- runif(1, min=0.6, max=0.9) # runif(1, min=0.1, max=0.5)
      oddscontrol <- pcontrol/(1-pcontrol)
      oddstreat <- oddscontrol * curor
      ptreat <- oddstreat/(1+oddstreat)


      control<-rbinom(cur_indnum0, 1, pcontrol)
      treat<-rbinom(cur_indnum1, 1, ptreat)

      datay<-c(control,treat)
      datax<-c(rep(0,cur_indnum0),rep(1,cur_indnum1))
      est_p<-summary(glm(datay~datax,
                         family="binomial"))$coefficient[2, c(1,2,4)]
      # if (rbinom(1,1,wipi(est_p[3])) == 1){  # || est_p[3] <= 0.05
      #   betalist<-c(betalist,est_p[1])
      #   sdlist<-c(sdlist,est_p[2])
      #   break
      # } else {
      #   if (est_p[3] > 0.1) { pubbias[jj, j] <- 1 } #  && est_p[3] <= 0.5
      # }

      if (rbinom(1,1,wipi(est_p[3], c)) == 1){  # || est_p[3] <= 0.05
        betalist<-c(betalist,est_p[1])
        sdlist<-c(sdlist,est_p[2])
        break
      } else {
        pubbias[jj, j] <- 1
      }
    }
  }

  betafinal <- rbind(betafinal,betalist)
  sdfinal <- rbind(sdfinal,sdlist)
}


hat_beta <- as.matrix(betafinal)
hat_sigma_sq <- as.matrix(sdfinal)^2
pubbiass <- rowSums(pubbias)
table(pubbiass)

save(hat_beta,
     file = paste0("../output/hat_beta_c_", c, ".rda"))
save(hat_sigma_sq,
     file = paste0("../output/hat_sigma_sq_c_", c, ".rda"))
save(pubbiass,
     file = paste0("../output/pubbiass_c_", c, ".rda"))

# pvec <- seq(0, 0.05, length = 1000)
# k_vec <- sapply(pvec, inverse_P_mis)
pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
# load("D:/R/Replicability/Peng/batch effect/data/k_vec.rda")
N <- 10000
r <- 0.05 # burn-in rate

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)
sim_results <- foreach(x = 1:num_simulations, .combine = combine_results, .packages = "DiscRep") %dopar% {
  metropolis_hastings(N, r, m, hat_beta[x,], hat_sigma_sq[x,], test = "Egger", k_vec = k_vec) #
}
stopCluster(cl)

# get p-value and bar_beta
p_values <- sim_results$p_value
bar_betas <- sim_results$bar_beta
sen <- sum(p_values <= 0.05)/length(p_values)


save(p_values,
     file = paste0("../output/pvalue_c_", c, ".rda"))

ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(0.00, 1.00)) +
  theme_bw() +
  theme(plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion")

# hist(bar_betas)
sen
sim_results[["acception"]]
max(p_values)
mean(p_values)
t2 <- Sys.time()
(t2-t1)





# # Egger
# rm(list = setdiff(ls(), "c"))
# # set.seed(23)
# library(ggplot2)
# library(foreach)
# library(doParallel)
# library(metafor)
# load(paste0("../output/hat_beta_c_", c, ".rda"))
# load(paste0("../output/hat_sigma_sq_c_", c, ".rda"))
#
# num_simulations <- nrow(hat_beta)
# EggerTest <- function(hat_beta, hat_sigma_sq) {
#   data <- data.frame(y = hat_beta, s2 = hat_sigma_sq)
#   res <- rma(yi = y, sei = sqrt(s2), data = data, method = "FE") #
#   egger_test <- regtest(res, model = "lm") #
#   return(egger_test$pval)
# }
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# p_values <- foreach(x = 1:num_simulations, .combine = c, .packages=c("metafor")) %dopar% {
#   EggerTest(hat_beta[x,], hat_sigma_sq[x,])
# }
# stopCluster(cl)
#
# sum(p_values <= 0.1)/length(p_values)
#
# save(p_values,
#      file = paste0("../output/pvalue_cEgger_", c, ".rda"))
# ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
#   geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
#   scale_x_continuous(limits = c(0.00, 1.00)) +
#   theme_bw() +
#   theme(plot.title = element_text(size=12, hjust = 0.5)) +
#   xlab("Posterior p-values under distinguishability criterion")
#
#
#
#
#
#
#
# # Begg
# rm(list = setdiff(ls(), "c"))
# # set.seed(23)
# library(ggplot2)
# library(foreach)
# library(doParallel)
# library(metafor)
# load(paste0("../output/hat_beta_c_", c, ".rda"))
# load(paste0("../output/hat_sigma_sq_c_", c, ".rda"))
#
#
# num_simulations <- nrow(hat_beta)
# BeggTest <- function(hat_beta, hat_sigma_sq) {
#   data <- data.frame(y = hat_beta, s2 = hat_sigma_sq)
#   res <- rma(yi = y, sei = sqrt(s2), data = data, method = "FE")
#   begg_test <- ranktest(res)
#   return(begg_test$pval)
# }
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# p_values <- foreach(x = 1:num_simulations, .combine = c, .packages=c("metafor")) %dopar% {
#   BeggTest(hat_beta[x,], hat_sigma_sq[x,])
# }
# stopCluster(cl)
# sum(p_values <= 0.1)/length(p_values)
#
#
#
# save(p_values,
#      file = paste0("../output/pvalue_cBegg_", c, ".rda"))
# ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
#   geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
#   scale_x_continuous(limits = c(0.00, 1.00)) +
#   theme_bw() +
#   theme(plot.title = element_text(size=12, hjust = 0.5)) +
#   xlab("Posterior p-values under distinguishability criterion")










# # Q-based posterior-PRP
rm(list = setdiff(ls(), "c"))
library(DiscRep)

load(paste0("../output/hat_beta_c_", c, ".rda"))
load(paste0("../output/hat_sigma_sq_c_", c, ".rda"))

set.seed(234)
num_simulations <- nrow(hat_beta)
studynum <- c(10, 6, 4)
m <- sum(studynum)

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
sen <- sum(p_values <= 0.05)/length(p_values)


save(p_values,
     file = paste0("../output/pvalue_cQ_", c, ".rda"))
ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(0.00, 1.00)) +
  theme_bw() +
  theme(plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion")
