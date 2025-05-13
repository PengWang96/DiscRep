rm(list = ls())
t1 <- Sys.time()
library(metafor)
library(DiscRep)
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


# # Stead et al. (2012)
data("dat.slf")
data <- dat.slf

# # Hrobjartsson and Gøtzsche (2010)
# data("dat.ha")
# data <- dat.ha

# # # Liu CJ, Latham NK (2009)
# data("dat.lcj")
# data <- dat.lcj



res <- rma(yi = y, sei = sqrt(s2), data = data) # , method = "FE"
res
par(mar = c(4.5, 5.5, 2, 2) + 0.1)
par(cex.axis = 1.5, cex.lab = 1.9) # cex.axis 控制坐标轴刻度的字体大小，cex.lab 控制坐标轴标签的字体大小
funnel(res)
# ggsave("./plot/funnel_slf.pdf", plot)


egger_test <- regtest(res, model = "lm", predictor = "sei")
print(egger_test)


begg_test <- ranktest(res)
print(begg_test)




# pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
# k_vec <- sapply(pvec, inverse_P_mis)
pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
m <- nrow(data)
hat_beta <- data$y
hat_sigma_sq <- data$s2
sim_results <- metropolis_hastings(10000, 0.05, m,
                                   hat_beta, hat_sigma_sq,
                                   test = "Egger", k_vec = k_vec) # Egger
(sim_results$p_value)



# Q test
frequency_pvalue(hat_beta, hat_sigma_sq)







library(foreach)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
p_values <- numeric(100)
p_values <- foreach(i=1:100, .combine='c', .packages = "DiscRep") %dopar% {
  sim_results <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq,
                                     test = "Egger", k_vec = k_vec) #
  return(sim_results$p_value)
}
stopCluster(cl)

estimate <- mean(p_values)
conf_interval <- quantile(p_values, c(0.025, 0.975))
print(paste("Estimate: ", estimate))
print(paste("95% Confidence Interval: ", conf_interval[1], ", ", conf_interval[2]))

t2 <- Sys.time()
(t2-t1)
