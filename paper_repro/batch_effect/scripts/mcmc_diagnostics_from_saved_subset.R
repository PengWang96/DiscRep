rm(list = ls())
library(ggplot2)
library(patchwork)
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
source("../../plot_style.R")

style <- plot_fullwidth_style(
  base_size = 8,
  width = 7,
  height = 4.5,
  grid = TRUE
)

t1 <- Sys.time()
set.seed(234)


bb_sd_vec <- c(0, 0.4, 0.8)
num_simulations <- 2000
n_pick_per_bb <- 4
use_common_sim_ids <- TRUE
r <- 0.05
max_lag <- 50
trace_thin <- 15

out_dir <- "../output/mcmc_diagnostics_saved_subset"
plot_dir <- "../plot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

compute_ess <- function(x, max_lag = 1000) {
  n <- length(x)
  if (n < 3) {
    return(NA_real_)
  }
  lag_max <- min(max_lag, n - 1)
  ac <- as.numeric(stats::acf(x, lag.max = lag_max, plot = FALSE, demean = TRUE)$acf)[-1]
  if (!length(ac)) {
    return(n)
  }
  first_negative <- which(ac < 0)
  if (length(first_negative)) {
    ac <- ac[seq_len(first_negative[1] - 1)]
  }
  tau <- 1 + 2 * sum(ac)
  if (!is.finite(tau) || tau <= 0) {
    return(n)
  }
  max(1, min(n, n / tau))
}

get_autocorr <- function(x, lag_value) {
  lag_max <- min(lag_value, length(x) - 1)
  if (lag_max < 1) {
    return(NA_real_)
  }
  ac <- stats::acf(x, lag.max = lag_max, plot = FALSE, demean = TRUE)$acf
  as.numeric(ac[lag_max + 1])
}

estimate_acceptance <- function(chain_draws) {
  if (length(chain_draws) < 2) {
    return(NA_real_)
  }
  mean(diff(chain_draws) != 0)
}

make_tag <- function(x) {
  gsub("\\.", "_", as.character(x))
}

diag_rows <- list()
trace_rows <- list()
acf_rows <- list()
picked_rows <- list()
all_chain_rows <- list()
diag_id <- 0
trace_id <- 0
acf_id <- 0
pick_id <- 0
all_chain_id <- 0

common_picked_idx <- NULL
if (use_common_sim_ids) {
  common_picked_idx <- sort(sample.int(num_simulations, n_pick_per_bb))
}

for (bb_sd in bb_sd_vec) {
  sim_path <- paste0("../output/sim_results_eta_", bb_sd, ".rds")
  if (!file.exists(sim_path)) {
    stop("Missing file: ", sim_path)
  }
  sim_results <- readRDS(sim_path)

  if (!is.matrix(sim_results$bar_beta)) {
    stop("Expected sim_results$bar_beta to be a matrix in file: ", sim_path)
  }
  if (nrow(sim_results$bar_beta) != num_simulations) {
    stop("Expected ", num_simulations, " rows in bar_beta for file: ", sim_path)
  }

  picked_idx <- if (use_common_sim_ids) common_picked_idx else sort(sample.int(num_simulations, n_pick_per_bb))
  pick_id <- pick_id + 1
  picked_rows[[pick_id]] <- data.frame(
    bb_sd = bb_sd,
    selected_sim_id = picked_idx
  )

  n_iter_all <- ncol(sim_results$bar_beta)
  burn_in_all <- floor(n_iter_all * r)
  if (burn_in_all >= (n_iter_all - 1)) {
    stop("Burn-in leaves too few draws for bb_sd = ", bb_sd)
  }
  post_idx_all <- (burn_in_all + 1):n_iter_all
  post_draws_all <- sim_results$bar_beta[, post_idx_all, drop = FALSE]
  all_posterior_mean <- rowMeans(post_draws_all)
  all_ess <- apply(
    post_draws_all,
    1,
    function(z) compute_ess(as.numeric(z))
  )
  all_chain_acceptance <- apply(
    sim_results$bar_beta,
    1,
    function(z) estimate_acceptance(as.numeric(z))
  )

  all_chain_id <- all_chain_id + 1
  all_chain_rows[[all_chain_id]] <- data.frame(
    bb_sd = bb_sd,
    posterior_mean = mean(all_posterior_mean, na.rm = TRUE),
    posterior_mean_sd_all2000 = stats::sd(all_posterior_mean, na.rm = TRUE),
    ess = mean(all_ess, na.rm = TRUE),
    chain_acceptance = mean(all_chain_acceptance, na.rm = TRUE),
    ess_sd_all2000 = stats::sd(all_ess, na.rm = TRUE),
    chain_acceptance_sd_all2000 = stats::sd(all_chain_acceptance, na.rm = TRUE),
    n_simulations_for_ess_acceptance = length(all_ess)
  )

  for (sim_id in picked_idx) {
    chain_full <- as.numeric(sim_results$bar_beta[sim_id, ])
    n_iter <- length(chain_full)
    burn_in <- floor(n_iter * r)
    post_draws <- chain_full[(burn_in + 1):n_iter]
    post_iter <- (burn_in + 1):n_iter

    diag_id <- diag_id + 1
    diag_rows[[diag_id]] <- data.frame(
      bb_sd = bb_sd,
      sim_id = sim_id,
      n_iter = n_iter,
      burn_in = burn_in,
      posterior_mean = mean(post_draws),
      posterior_sd = stats::sd(post_draws),
      ess = compute_ess(post_draws),
      acf_lag1 = get_autocorr(post_draws, 1),
      acf_lag5 = get_autocorr(post_draws, 5),
      acf_lag10 = get_autocorr(post_draws, 10),
      chain_acceptance = estimate_acceptance(chain_full),
      p_value = sim_results$p_value[sim_id]
    )

    keep <- seq.int(1, length(post_draws), by = trace_thin)
    trace_id <- trace_id + 1
    trace_rows[[trace_id]] <- data.frame(
      bb_sd = bb_sd,
      sim_id = sim_id,
      iter = post_iter[keep],
      bar_beta = post_draws[keep]
    )

    ac <- stats::acf(post_draws, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf
    acf_id <- acf_id + 1
    acf_rows[[acf_id]] <- data.frame(
      bb_sd = bb_sd,
      sim_id = sim_id,
      lag = 0:max_lag,
      acf = as.numeric(ac)
    )
  }
  cat("Completed diagnostics from saved chains for bb_sd =", bb_sd, "\n")
}

picked_df <- do.call(rbind, picked_rows)
diag_df <- do.call(rbind, diag_rows)
trace_df <- do.call(rbind, trace_rows)
acf_df <- do.call(rbind, acf_rows)
all_chain_df <- do.call(rbind, all_chain_rows)

diag_summary <- all_chain_df

write.csv(picked_df, file.path(out_dir, "selected_chain_indices.csv"), row.names = FALSE)
write.csv(diag_df, file.path(out_dir, "selected_chain_diagnostics.csv"), row.names = FALSE)
write.csv(diag_summary, file.path(out_dir, "selected_chain_diagnostics_summary.csv"), row.names = FALSE)

eta_expr_levels <- paste0("eta == ", bb_sd_vec)
eta_lookup <- data.frame(
  bb_sd = bb_sd_vec,
  eta_expr = eta_expr_levels
)
trace_plot_df <- merge(trace_df, eta_lookup, by = "bb_sd", sort = FALSE)
acf_plot_df <- merge(acf_df, eta_lookup, by = "bb_sd", sort = FALSE)
trace_plot_df$sim_id <- factor(trace_plot_df$sim_id, levels = sort(unique(trace_plot_df$sim_id)))
acf_plot_df$sim_id <- factor(acf_plot_df$sim_id, levels = sort(unique(acf_plot_df$sim_id)))
trace_plot_df$eta_expr <- factor(trace_plot_df$eta_expr, levels = eta_expr_levels)
acf_plot_df$eta_expr <- factor(acf_plot_df$eta_expr, levels = eta_expr_levels)

p_trace <- ggplot(trace_plot_df, aes(x = iter, y = bar_beta, color = sim_id, group = sim_id)) +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  facet_wrap(
    ~eta_expr,
    nrow = 1,
    labeller = labeller(eta_expr = label_parsed)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(
    limits = c(501, 10000),
    breaks = c(501, 2500, 5000, 7500, 10000)
  ) +
  labs(
    x = "Iteration",
    y = expression(bar(beta))
  ) +
  style$theme +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "plain"),
    strip.text = element_text(face = "plain", size = 8.5),
    plot.margin = margin(5, 6, 5, 6)
  )

p_acf <- ggplot(acf_plot_df, aes(x = lag, y = acf, color = sim_id, group = sim_id)) +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.35) +
  facet_wrap(
    ~eta_expr,
    nrow = 1,
    labeller = labeller(eta_expr = label_parsed)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(
    limits = c(0, 50),
    breaks = c(0, 10, 20, 30, 40, 50)
  ) +
  labs(
    x = "Lag",
    y = "ACF"
  ) +
  style$theme +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "plain"),
    strip.text = element_text(face = "plain", size = 8.5),
    plot.margin = margin(5, 6, 5, 6)
  )

p_facet <- p_trace / p_acf + plot_layout(heights = c(1, 1))

style$save_pdf(
  p_facet,
  file.path(plot_dir, "trace_acf_facet_by_eta"),
  w = 7,
  h = 4.5
)

cat("\nSaved diagnostics to:", normalizePath(out_dir), "\n")
cat("Saved plots to:", normalizePath(plot_dir), "\n")
t2 <- Sys.time()
cat("\nTotal elapsed time:", format(t2 - t1), "\n")
