##############################################################################
######## Histogram: p-values under different selection levels
##############################################################################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(rstudioapi)

if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

source("../../plot_style.R")

style <- plot_fullwidth_style(
  base_size = 6,
  width = 7,
  height = 4.5,
  grid = TRUE
)

c_levels <- c(0, 10, 20)
censor_levels <- c(
  "No selection (c = 0)",
  "Modest selection (c = 10)",
  "Strong selection (c = 20)"
)
test_levels <- c(
  "Modified Egger statistic",
  "Modified Q statistic"
)

# Load p-values for Modified Egger and Modified Q statistics.
egger_values <- lapply(c_levels, function(c_val) {
  readRDS(paste0("../output/pvalue_c_", c_val, ".rds"))
})
q_values <- lapply(c_levels, function(c_val) {
  readRDS(paste0("../output/pvalue_cQ_", c_val, ".rds"))
})

# Build long-format data for 2x3 facets (rows = tests, columns = censor levels).
combined_data <- bind_rows(
  bind_rows(lapply(seq_along(c_levels), function(i) {
    data.frame(
      Censor = censor_levels[i],
      Test = "Modified Egger statistic",
      P_Value = egger_values[[i]]
    )
  })),
  bind_rows(lapply(seq_along(c_levels), function(i) {
    data.frame(
      Censor = censor_levels[i],
      Test = "Modified Q statistic",
      P_Value = q_values[[i]]
    )
  }))
)

combined_data$Censor <- factor(combined_data$Censor, levels = censor_levels)
combined_data$Test <- factor(combined_data$Test, levels = test_levels)

# Show mean reference only in the c = 0 panels.
means_data <- combined_data %>%
  filter(Censor == "No selection (c = 0)") %>%
  group_by(Censor, Test) %>%
  summarize(mean_p = mean(P_Value), .groups = "drop")
print(means_data)

plot <- ggplot(combined_data, aes(x = P_Value)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  facet_grid(Test ~ Censor) +
  geom_vline(
    data = means_data,
    aes(xintercept = mean_p),
    color = "blue",
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = means_data,
    aes(x = mean_p, y = Inf, label = sprintf("Mean: %.3f", mean_p)),
    color = "blue",
    family = style$base_family,
    vjust = 2,
    hjust = -0.1,
    size = style$geom_text_size_mm(6),
    inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  labs(
    x = "Posterior-PRPs",
    y = "Count"
  ) +
  style$theme +
  theme(
    axis.title = element_text(face = "plain"),
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    strip.text = element_text(face = "plain", size = 6.5),
    plot.margin = margin(5, 6, 5, 6)
  )
plot
style$save_pdf(plot, "../plot/pub_bias_pvalue", w = 7, h = 4.5)














##############################################################################
######## Sensitivity: compare methods under different selection levels
##############################################################################
rm(list = ls())
library(ggplot2)
data <- data.frame(
  c = c(0, 2.5, 5, 7.5, 10),
  `Posterior-PRP` = c(0.095, 0.165, 0.185, 0.255, 0.29),
  Egger = c(0.08, 0.13, 0.13, 0.15, 0.215),
  Begg = c(0.085, 0.085, 0.115, 0.2, 0.21)
)
data_long <- reshape2::melt(data, id.vars = "c")
ggplot(data = data_long, aes(x = c, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(x = expression(eta), y = "Sensitivity") +
  theme_bw() +
  scale_y_continuous(limits = c(0.05, 0.4)) +
  scale_color_discrete(labels = c("Posterior-PRP", "Egger Test", "Begg Test")) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.direction = "vertical"
  )


rm(list = ls())
library(ggplot2)
data <- data.frame(
  c = c(0, 2.5, 5, 7.5, 10),
  `Posterior-PRP` = c(0.121, 0.1275, 0.185, 0.249, 0.27),
  Egger = c(0.1285, 0.139, 0.139, 0.162, 0.162),
  Begg = c(0.0795, 0.0875, 0.1295, 0.1765, 0.2095)
)
data_long <- reshape2::melt(data, id.vars = "c")
ggplot(data = data_long, aes(x = c, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(x = "c", y = "Sensitivity") +
  theme_bw() +
  scale_y_continuous(limits = c(0.05, 0.3)) +
  scale_color_discrete(labels = c("Posterior-PRP", "Egger Test", "Begg Test")) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.direction = "vertical"
  )
