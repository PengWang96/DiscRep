######## histogram: 不同 eta 对应的 p value
rm(list = ls())
library(ggplot2)
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
bb_sd <- 0
load(file = paste0("../output/pvalue_eta_", bb_sd, ".rda"))
mean_p <- mean(p_values)
pval <- data.frame(p_value = p_values)

sum(p_values <= 0.05)/length(p_values)

plot <- ggplot(pval, aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  # scale_x_continuous(limits = c(-0.01, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 15), # Increase x-axis title font size
    axis.title.y = element_text(size = 15)  # Increase y-axis title font size
  ) +
  geom_vline(data = pval, aes(xintercept = mean_p), color = "blue", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = mean_p, y = 200, label = paste0("Mean: ", round(mean_p, 3)),
           color = "blue", angle = 0, vjust = -2, hjust = -0.1, size = 5) +
  xlab("Posterior-PRPs under Distinguishability Criterion") +
  ylab("Count") +
  ylim(0, 1500)
plot
ggsave("../plot/batch_eta_0.pdf", plot)
dev.off()



rm(list = ls())
library(ggplot2)
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
bb_sd <- 0.4 #0.8 #
load(file = paste0("../output/pvalue_eta_", bb_sd, ".rda"))
pval <- data.frame(p_value = p_values)
2 * mean(p_values)
sum(p_values <= 0.05)/length(p_values)

plot <- ggplot(pval, aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  # scale_x_continuous(limits = c(-0.02, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 15), # Increase x-axis title font size
    axis.title.y = element_text(size = 15)  # Increase y-axis title font size
  ) +
  xlab("Posterior-PRPs under Distinguishability Criterion") +
  ylab("Count") +
  ylim(0, 1500)
plot
ggsave(paste0("../plot/batch_eta_", bb_sd, ".pdf"), plot)
dev.off()



######## Sensitivity: 不同 eta 对应的 Sensitivity
rm(list = ls())
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
library(ggplot2)
data <- data.frame(eta = c(0, 0.4, 0.8, 1.2),
                   Sensitivity = c(0.0475, 0.306, 0.7485, 0.909))
data_long <- reshape2::melt(data, id.vars = "eta")
plot <- ggplot(data=data_long, aes(x=eta, y=value)) +
  geom_line(linewidth=1, color = "#6186ad") +
  geom_point(size=3, color = "#6186ad") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Batch Contamination Level", y = "Proportion of Flagged Datasets") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 15), # Increase x-axis title font size
    axis.title.y = element_text(size = 15)  # Increase y-axis title font size
  ) +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
                     limits = c(0.00, 1))
plot
ggsave("../plot/Sensitivity.pdf", plot)
dev.off()




######## histogram: fix effect 的 p value. 基于Frequentist, 理论上的，还有算法1得到的
rm(list = ls())
library(ggplot2)
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
bb_sd <- 0
# load(file = paste0("../output/Frequentist_pvalue_eta_", bb_sd, ".rda"))
load(file = "../output/pvalue_eta_0_k_0.rda")

mean_p <- mean(p_values)

pval <- data.frame(p_value = p_values)
plot <- ggplot(pval, aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(0.00, 1.00)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 17), # Increase x-axis title font size
    axis.title.y = element_text(size = 17)  # Increase y-axis title font size
  ) +
  geom_vline(data = pval, aes(xintercept = mean_p), color = "blue", linetype = "dashed", size = 0.6) +
  annotate("text", x = mean_p, y = 200, label = paste0("Mean: ", round(mean_p, 3)),
           color = "blue", angle = 0, vjust = 1.5, hjust = -0.1, size = 7) +
  # xlab("Frequentist p-values") +
  xlab("Posterior-PRPs") +
  ylab("Count") +
  ylim(0, 200)
plot
# ggsave("../plot/fixed_effect_frequentist_pvalues.pdf", plot)
ggsave("../plot/fixed_effect_algorithm_pvalues.pdf", plot)
dev.off()


































######## Sensitivity: 不同 eta 对应的 不同方法的Sensitivity
# bbar <- 1
rm(list = ls())
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
library(ggplot2)
data <- data.frame(eta = c(0, 0.4, 0.8, 1.2),
                   `Posterior-PRP` = c(0.615, 0.7, 0.825, 0.92),
                   `Q Test` = c(0.5, 0.635, 0.785, 0.895))
data_long <- reshape2::melt(data, id.vars = "eta")
ggplot(data=data_long, aes(x=eta, y=value, color = variable)) +
  geom_line(linewidth=1) +
  geom_point(size=3) +
  labs(x = expression(eta), y = "Sensitivity") +
  theme_bw() +
  scale_y_continuous(limits = c(0.4, 0.95)) +
  scale_color_discrete(labels = c("Posterior-PRP", "Q Test")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical")




######## 检查MCMC收敛性的几种图
# # 画出迭代次数和 bar_beta 的折线图
for (i in 1:10) {
  p <- ggplot(data.frame(iteration = 9000:N,
                    bar_beta = bar_betas[i,9000:N]),
         aes(x = iteration, y = bar_beta)) +
    geom_line() +
    labs(title = 'MH MCMC Convergence',
         x = 'Iteration',
         y = expression('Sampled' ~ bar(beta)))

  print(p)
  Sys.sleep(1)  # Pause for visualization
}

hist(bar_betas)

for (i in 1:5) {
  acf(bar_betas[i, (N*r):N])
  Sys.sleep(1)  # Pause for visualization
}



