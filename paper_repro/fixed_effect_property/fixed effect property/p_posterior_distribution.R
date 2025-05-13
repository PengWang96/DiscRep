rm(list = ls())
library(ggplot2)
library(dplyr)

compute_p_posterior <- function(c0, m) {
  integrand <- function(s) {
    one_minus_cdf <- 1 - pchisq(s + c0, df = m)
    pdf <- dchisq(s, df = 1)
    return(one_minus_cdf * pdf)
  }

  result <- integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-8)

  return(result$value)
}


m_values <- c(2, 4, 6, 8, 10, 12)
num_simulations <- 100000
set.seed(123)

simulation_results <- lapply(m_values, function(m) {
  true_beta <- 0
  sigma_sq <- rep(0.01, m)
  w_j <- 1 / sigma_sq

  # \hat{\beta}_j ~ N(\bar{\beta}, \hat{\sigma}_j^2)
  beta_hat_matrix <- matrix(rnorm(m * num_simulations, mean = true_beta, sd = sqrt(sigma_sq)),
                            nrow = num_simulations, ncol = m)

  # \mu_{\bar{\beta}} = sum(w_j * \hat{\beta}_j) / sum(w_j)
  mu_bar_beta <- rowSums(beta_hat_matrix * matrix(w_j, nrow = num_simulations, ncol = m, byrow = TRUE)) / sum(w_j)

  # c0 = sum(w_j * (\hat{\beta}_j - mu_bar_beta)^2)
  c0 <- rowSums(w_j * (beta_hat_matrix - mu_bar_beta)^2)

  # p_posterior
  p_posterior <- sapply(c0, compute_p_posterior, m = m)

  data.frame(m = m, c0 = c0, p_posterior = p_posterior)
}) %>% bind_rows()


summary(simulation_results$p_posterior)
sum(simulation_results$p_posterior <= 0)
sum(simulation_results$p_posterior > 1)


# plot p_posterior vs. c0
ggplot(simulation_results, aes(x = c0, y = p_posterior, color = as.factor(m))) +
  geom_line(size = 1) +
  labs(
    title = expression(paste(p[posterior-PRP], " vs. ", c[0])),
    x = expression(c[0]),
    y = expression(p[posterior-PRP]),
    color = "m"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )



summary_stats <- simulation_results %>%
  group_by(m) %>%
  summarise(
    mean_p = mean(p_posterior)
  )

ggplot(simulation_results, aes(x = p_posterior, fill = as.factor(m))) +
  geom_histogram(bins = 100, alpha = 0.8, position = "identity", color = "white") +
  facet_wrap(~ m, nrow = 3, scales = "free_y") +
  labs(
    title = NULL,
    x = expression(p[posterior-PRP]),
    y = "Count",
    fill = "m"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  ) +
  geom_vline(data = summary_stats, aes(xintercept = mean_p), color = "blue", linetype = "dashed", size = 0.6) +
  geom_text(
    data = summary_stats,
    aes(x = mean_p, y = Inf, label = paste0("Mean: ", round(mean_p, 3))),
    color = "blue",
    angle = 0,
    vjust = 2,
    hjust = 2,
    size = 5
  ) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ m, nrow = 3, scales = "free_y") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.ticks.x = element_line()
  )



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





