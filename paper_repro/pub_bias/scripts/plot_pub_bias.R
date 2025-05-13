######## histogram: 不同 c 对应的 p value
rm(list = ls())
library(ggplot2)
library(dplyr)
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

load(file = paste0("../output/pvalue_c_", 0, ".rda"))
p_values0 <- p_values
load(file = paste0("../output/pvalue_c_", 10, ".rda"))
p_values10 <- p_values; 2 * mean(p_values10)
load(file = paste0("../output/pvalue_c_", 20, ".rda"))
p_values20 <- p_values; 2 * mean(p_values20)
load(file = paste0("../output/pvalue_cQ_", 0, ".rda"))
p_valuesQ0 <- p_values
load(file = paste0("../output/pvalue_cQ_", 10, ".rda"))
p_valuesQ10 <- p_values
load(file = paste0("../output/pvalue_cQ_", 20, ".rda"))
p_valuesQ20 <- p_values

# 创建组合数据框
combined_data <- data.frame(
  Censor = rep(c("No Selection (c = 0)", "Modest Selection (c = 10)", "Strong Selection (c = 20)"),
               each = length(p_valuesQ0),
               times = 2),
  Test = rep(c("Modified Q Statistic", "Modified Egger Statistic"), each = length(p_valuesQ0) * 3),
  P_Value = c(p_valuesQ0, p_valuesQ10, p_valuesQ20, p_values0, p_values10, p_values20)
)

combined_data$Censor <- factor(combined_data$Censor,
                               levels = c("No Selection (c = 0)",
                                          "Modest Selection (c = 10)",
                                          "Strong Selection (c = 20)"))
head(combined_data)

# 计算均值（仅针对 Censor = "No censoring (c = 0)"）
means_data <- combined_data %>%
  filter(Censor == "No Selection (c = 0)") %>%
  group_by(Censor, Test) %>%
  summarize(mean_p = mean(P_Value))
print(means_data)


plot <- ggplot(combined_data, aes(x = P_Value)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  facet_grid(Censor ~ Test) +
  labs(
    x = "Posterior-PRPs under Distinguishability Criterion",
    y = "Count"
  ) +
  theme_bw() +
  # 添加垂直线，仅在 Censor = "No censoring (c = 0)" 的子图中
  geom_vline(data = means_data, aes(xintercept = mean_p), color = "blue",
             linetype = "dashed", linewidth = 0.6) +
  # 添加均值标签，仅在 Censor = "No censoring (c = 0)" 的子图中
  geom_text(
    data = means_data,
    aes(x = mean_p, y = Inf, label = sprintf("Mean: %.3f", mean_p)),
    color = "blue",
    vjust = 1.5,
    hjust = -0.1,
    size = 2.8,
    inherit.aes = FALSE
  )  +
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  theme(
    axis.title.x = element_text(size = 8), # Increase x-axis title font size
    axis.title.y = element_text(size = 8),  # Increase y-axis title font size
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
plot
ggsave("../plot/pub_bias_pvalue.pdf", plot, width = 6, height = 6)
dev.off()














######## Sensitivity: 不同 eta 对应的 不同方法的Sensitivity
rm(list = ls())
library(ggplot2)
data <- data.frame(c = c(0, 2.5, 5, 7.5, 10),
                   `Posterior-PRP` = c(0.095, 0.165, 0.185, 0.255, 0.29),
                   Egger = c(0.08, 0.13, 0.13, 0.15, 0.215),
                   Begg = c(0.085, 0.085, 0.115, 0.2, 0.21))
data_long <- reshape2::melt(data, id.vars = "c")
ggplot(data=data_long, aes(x=c, y=value, color = variable)) +
  geom_line(linewidth=1) +
  geom_point(size=3) +
  labs(x = expression(eta), y = "Sensitivity") +
  theme_bw() +
  scale_y_continuous(limits = c(0.05, 0.4)) +
  scale_color_discrete(labels = c("Posterior-PRP", "Egger Test", "Begg Test")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical")



rm(list = ls())
library(ggplot2)
data <- data.frame(c = c(0, 2.5, 5, 7.5, 10),
                   `Posterior-PRP` = c(0.121, 0.1275, 0.185, 0.249, 0.27),
                   Egger = c(0.1285, 0.139, 0.139, 0.162, 0.162),
                   Begg = c(0.0795, 0.0875, 0.1295, 0.1765, 0.2095))
data_long <- reshape2::melt(data, id.vars = "c")
ggplot(data=data_long, aes(x=c, y=value, color = variable)) +
  geom_line(linewidth=1) +
  geom_point(size=3) +
  labs(x = "c", y = "Sensitivity") +
  theme_bw() +
  scale_y_continuous(limits = c(0.05, 0.3)) +
  scale_color_discrete(labels = c("Posterior-PRP", "Egger Test", "Begg Test")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical")

