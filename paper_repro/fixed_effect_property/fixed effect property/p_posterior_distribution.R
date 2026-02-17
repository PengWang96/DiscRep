rm(list = ls())
library(ggplot2)
library(dplyr)
library(rstudioapi)
t1 <- Sys.time()
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

source("../../plot_style.R")

# Keep typography aligned with batch_effect/scripts/plot_batch.R.
style <- plot_fullwidth_style(
  base_size = 10,
  width = 7,
  height = 5.6,
  grid = TRUE
)
p_vs_c0_font_pt <- style$base_size + 1.5
prp_hist_font_pt <- style$base_size - 1.2

plot_dir <- "./plot"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

compute_p_posterior <- function(c0, m) {
  integrand <- function(s) {
    one_minus_cdf <- 1 - pchisq(s + c0, df = m)
    pdf <- dchisq(s, df = 1)
    one_minus_cdf * pdf
  }

  result <- integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-8)
  result$value
}

m_values <- c(2, 4, 6, 8, 10, 12)
num_simulations <- 100000
set.seed(123)

simulation_results <- lapply(m_values, function(m) {
  true_beta <- 0
  sigma_sq <- rep(0.01, m)
  w_j <- 1 / sigma_sq

  # \hat{\beta}_j ~ N(\bar{\beta}, \hat{\sigma}_j^2)
  beta_hat_matrix <- matrix(
    rnorm(m * num_simulations, mean = true_beta, sd = sqrt(sigma_sq)),
    nrow = num_simulations,
    ncol = m
  )

  # \mu_{\bar{\beta}} = sum(w_j * \hat{\beta}_j) / sum(w_j)
  mu_bar_beta <- rowSums(
    beta_hat_matrix * matrix(w_j, nrow = num_simulations, ncol = m, byrow = TRUE)
  ) / sum(w_j)

  # c0 = sum(w_j * (\hat{\beta}_j - mu_bar_beta)^2)
  c0 <- rowSums(w_j * (beta_hat_matrix - mu_bar_beta)^2)

  # p_posterior
  p_posterior <- sapply(c0, compute_p_posterior, m = m)

  data.frame(m = m, c0 = c0, p_posterior = p_posterior)
}) %>%
  bind_rows()

simulation_results <- simulation_results %>%
  mutate(
    m = factor(m, levels = m_values),
    m_label = factor(paste0("m==", as.character(m)), levels = paste0("m==", m_values))
  )

summary(simulation_results$p_posterior)
sum(simulation_results$p_posterior <= 0)
sum(simulation_results$p_posterior > 1)

# plot p_posterior vs. c0
plot_p_vs_c0 <- ggplot(simulation_results, aes(x = c0, y = p_posterior, color = m)) +
  geom_line(linewidth = 0.9) +
  labs(
    title = NULL,
    x = expression(c[0]),
    y = "Posterior-PRPs",
    color = "m"
  ) +
  scale_color_brewer(palette = "Set2") +
  style$theme +
  theme(
    axis.title = element_text(size = p_vs_c0_font_pt * 1.1, face = "plain"),
    axis.title.x = element_text(size = p_vs_c0_font_pt * 1.1, family = style$math_family),
    axis.title.y = element_text(size = p_vs_c0_font_pt * 1.1),
    axis.text = element_text(size = p_vs_c0_font_pt * 0.9),
    legend.title = element_text(size = p_vs_c0_font_pt * 1.0, face = "bold"),
    legend.text = element_text(size = p_vs_c0_font_pt * 0.9),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.8),
      color = NA
    ),
    plot.margin = margin(5, 6, 5, 6)
  )

style$save_pdf(
  plot_p_vs_c0,
  file.path(plot_dir, "p_vs_c0"),
  w = 7,
  h = 6
)

summary_stats <- simulation_results %>%
  group_by(m, m_label) %>%
  summarise(mean_p = mean(p_posterior), .groups = "drop")

plot_prp_hist <- ggplot(simulation_results, aes(x = p_posterior, fill = m)) +
  geom_histogram(bins = 100, alpha = 0.8, position = "identity", color = "white") +
  geom_vline(
    data = summary_stats,
    aes(xintercept = mean_p),
    color = "blue",
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = summary_stats,
    aes(x = mean_p, y = Inf, label = paste0("Mean: ", round(mean_p, 3))),
    family = style$base_family,
    color = "blue",
    vjust = 2,
    hjust = 1.02,
    size = style$geom_text_size_mm(prp_hist_font_pt),
    inherit.aes = FALSE
  ) +
  facet_wrap(~m_label, nrow = 3, scales = "free_y", labeller = label_parsed) +
  labs(
    title = NULL,
    x = "Posterior-PRPs",
    y = "Count",
    fill = "m"
  ) +
  scale_fill_brewer(palette = "Set2") +
  style$theme +
  theme(
    axis.title = element_text(size = prp_hist_font_pt * 1.1, face = "plain"),
    axis.text = element_text(size = prp_hist_font_pt * 0.9),
    legend.position = "none",
    strip.text = element_text(
      family = style$math_family,
      size = prp_hist_font_pt,
      face = "plain"
    ),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.ticks.x = element_line()
  )

style$save_pdf(
  plot_prp_hist,
  file.path(plot_dir, "PRP_hist"),
  w = 7,
  h = 7
)

# Keep manuscript figure folder in sync when it exists.
bio_figure_dir <- normalizePath("../../../Bioinformatics/figure", winslash = "/", mustWork = FALSE)
if (dir.exists(bio_figure_dir)) {
  file.copy(
    file.path(plot_dir, "p_vs_c0.pdf"),
    file.path(bio_figure_dir, "p_vs_c0.pdf"),
    overwrite = TRUE
  )
  file.copy(
    file.path(plot_dir, "PRP_hist.pdf"),
    file.path(bio_figure_dir, "PRP_hist.pdf"),
    overwrite = TRUE
  )
}
t2 <- Sys.time()
print(paste("Total runtime:", round(difftime(t2, t1, units = "secs"), 2), "seconds"))


# # p_posterior hist
# ggplot(simulation_results, aes(x = p_posterior, fill = as.factor(m))) +
#   geom_histogram(bins = 100, alpha = 0.8, position = "identity", color = "white") +
#   facet_wrap(~ m, nrow = 3) +
#   labs(
#     title = NULL, # expression(paste("Histogram of ", p[posterior])),
#     x = expression(p[posterior]),
#     y = "Count",
#     fill = "m"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 15)
#   )
#
#
# ggplot(simulation_results, aes(x = p_posterior, fill = as.factor(m))) +
#   geom_histogram(bins = 100, alpha = 0.8, position = "identity", color = "white") +
#   facet_wrap(~ m, nrow = 3, scales = "free_y") +
#   labs(
#     title = NULL, # expression(paste("Histogram of ", p[posterior])),
#     x = expression(p[posterior]),
#     y = "Count",
#     fill = "m"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 15)
#   )

