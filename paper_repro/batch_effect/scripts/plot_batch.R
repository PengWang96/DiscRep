##############################################################################
######## Histogram + Sensitivity: combined 2x2 figure
##############################################################################
rm(list = ls())
library(ggplot2)
library(patchwork)
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

# For double-column insertion with LaTeX scaling 0.8:
# base_size 7.5 -> effective size ~9.6 pt in manuscript.
style <- plot_fullwidth_style(
  base_size = 7.5,
  width = 7,
  height = 6.2,
  grid = TRUE
)
x_axis_title_pt <- style$base_size * 1.1

bb_sd_vec <- c(0, 0.4, 0.8)
sim_files <- file.path("../output", paste0("sim_results_eta_", bb_sd_vec, ".rds"))
missing_files <- sim_files[!file.exists(sim_files)]
if (length(missing_files) > 0) {
  stop("Missing file(s): ", paste(missing_files, collapse = ", "))
}

pval_list <- lapply(sim_files, function(f) {
  sim_results <- readRDS(f)
  data.frame(p_value = sim_results$p_value)
})
names(pval_list) <- as.character(bb_sd_vec)

hist_bins <- 15
hist_breaks <- seq(0, 1, length.out = hist_bins + 1)
max_count <- max(vapply(pval_list, function(d) {
  max(hist(d$p_value, breaks = hist_breaks, plot = FALSE)$counts)
}, numeric(1)))
hist_ymax <- ceiling(max_count * 1.05 / 100) * 100
hist_ymin <- -0.08 * hist_ymax

hist_plots <- Map(function(bb_sd, pval) {
  mean_p <- mean(pval$p_value)
  p <- ggplot(pval, aes(x = p_value)) +
    geom_histogram(
      color = "white",
      fill = "#6186ad",
      bins = hist_bins,
      boundary = 0,
      alpha = 0.80
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0.05)
    ) +
    scale_y_continuous(
      limits = c(hist_ymin, hist_ymax),
      breaks = pretty(c(0, hist_ymax), n = 6),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "Posterior-PRPs", y = "Count") +
    style$theme +
    theme(
      axis.title = element_text(face = "plain"),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      plot.margin = margin(5, 6, 2, 6)
    )

  if (isTRUE(all.equal(bb_sd, 0))) {
    p <- p +
      geom_vline(
        xintercept = mean_p,
        color = "blue",
        linetype = "dashed",
        linewidth = 0.6
      ) +
      annotate(
        "text",
        x = mean_p,
        y = hist_ymax * 0.15,
        label = paste0("Mean: ", round(mean_p, 3)),
        family = style$base_family,
        color = "blue",
        hjust = -0.05,
        vjust = -0.4,
        size = style$geom_text_size_mm(8)
      )
  }
  p
}, bb_sd_vec, pval_list)

sensitivity_eta_vec <- c(0, 0.4, 0.8, 1.2)
sensitivity_files <- file.path("../output", paste0("sim_results_eta_", sensitivity_eta_vec, ".rds"))
missing_sensitivity_files <- sensitivity_files[!file.exists(sensitivity_files)]
if (length(missing_sensitivity_files) > 0) {
  stop("Missing sensitivity file(s): ", paste(missing_sensitivity_files, collapse = ", "))
}

sensitivity_values <- vapply(sensitivity_files, function(f) {
  sim_results <- readRDS(f)
  mean(sim_results$p_value <= 0.05)
}, numeric(1))

sensitivity_df <- data.frame(
  eta = sensitivity_eta_vec,
  Sensitivity = sensitivity_values
)

cat("Sensitivity by eta (p_value <= 0.05):\n")
print(sensitivity_df)

sensitivity_plot <- ggplot(sensitivity_df, aes(x = eta, y = Sensitivity)) +
  geom_line(linewidth = 1, color = "#6186ad") +
  geom_point(size = 2.6, color = "#6186ad") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_x_continuous(
    breaks = sensitivity_eta_vec,
    limits = c(min(sensitivity_eta_vec), max(sensitivity_eta_vec))
  ) +
  scale_y_continuous(
    breaks = c(0.05, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Batch contamination level",
    y = "Proportion of flagged datasets"
  ) +
  style$theme +
  theme(
    axis.title = element_text(face = "plain"),
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.margin = margin(5, 6, 2, 6)
  )

make_label_panel <- function(tag, eta = NULL, text = NULL) {
  p <- ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  if (!is.null(eta)) {
    p <- p +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0("(", tag, ") "),
        family = style$base_family,
        hjust = 1,
        size = style$geom_text_size_mm(x_axis_title_pt)
      ) +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0("eta==", eta),
        parse = TRUE,
        family = style$math_family,
        hjust = 0,
        size = style$geom_text_size_mm(x_axis_title_pt)
      )
  } else {
    p <- p +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0("(", tag, ") ", text),
        family = style$base_family,
        hjust = 0.5,
        size = style$geom_text_size_mm(x_axis_title_pt)
      )
  }
  p
}

panel_with_bottom_label <- function(main_plot, label_plot) {
  main_plot / label_plot + plot_layout(heights = c(1, 0.12))
}

panel_a <- panel_with_bottom_label(hist_plots[[1]], make_label_panel("a", eta = 0))
panel_b <- panel_with_bottom_label(hist_plots[[2]], make_label_panel("b", eta = 0.4))
panel_c <- panel_with_bottom_label(hist_plots[[3]], make_label_panel("c", eta = 0.8))
panel_d <- panel_with_bottom_label(
  sensitivity_plot,
  make_label_panel("d", text = "Proportion of flagged datasets")
)

combined_plot <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(2, 2, 2, 2))

style$save_pdf(
  combined_plot,
  "../plot/batch_eta_sensitivity_2x2",
  w = 7,
  h = 6.2
)




##############################################################################
######## Histogram: fixed-effect p-values (traditional Q test)
# and posterior-PRP from Algorithm 1
##############################################################################
rm(list = ls())
library(ggplot2)
library(patchwork)
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

# For LaTeX scaling = 1: keep native full-width sizing.
style <- plot_fullwidth_style(
  base_size = 10,
  width = 7,
  height = 3.4,
  grid = TRUE
)

dat <- readRDS("../output/pvalues_k_0.rds")
required_names <- c("Q_bb_sd_0", "PRP_bb_sd_0")
if (!all(required_names %in% names(dat))) {
  stop("Missing required fields in dat: ", paste(setdiff(required_names, names(dat)), collapse = ", "))
}

q_df <- data.frame(p_value = dat$Q_bb_sd_0)
prp_df <- data.frame(p_value = dat$PRP_bb_sd_0)

build_fixed_hist <- function(df, x_label) {
  mean_p <- mean(df$p_value)
  ggplot(df, aes(x = p_value)) +
    geom_histogram(
      color = "white",
      fill = "#6186ad",
      bins = 20,
      boundary = 0,
      alpha = 0.8
    ) +
    geom_vline(
      xintercept = mean_p,
      color = "blue",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = mean_p,
      y = Inf,
      label = paste0("Mean: ", round(mean_p, 3)),
      family = style$base_family,
      color = "blue",
      hjust = -0.05,
      vjust = 2,
      size = style$geom_text_size_mm(9)
    ) +
    scale_y_continuous(
      limits = c(0, 200),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = x_label, y = "Count") +
    style$theme +
    theme(
      axis.title = element_text(face = "plain"),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      plot.margin = margin(5, 6, 5, 6)
    )
}

plot_q <- build_fixed_hist(q_df, "Frequentist p-values")
plot_prp <- build_fixed_hist(prp_df, "Posterior-PRPs")

combined_fixed_hist <- (plot_q | plot_prp) +
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(2, 2, 2, 2))

style$save_pdf(
  combined_fixed_hist,
  "../plot/fixed_effect_Q_vs_PRP_hist_eta0",
  w = 7,
  h = 3.4
)




##############################################################################
# scatter plot for k = 0, fixed effect, compare posterior-PRP  
# and Q test based on batch effect simulation
##############################################################################
rm(list = ls())
library(ggplot2)
library(patchwork)
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

# For LaTeX scaling = 1: keep native full-width sizing.
style <- plot_fullwidth_style(
  base_size = 10,
  width = 7,
  height = 2.3,
  grid = TRUE
)
x_axis_title_pt <- style$base_size * 1.1

pvalue_results <- readRDS("../output/pvalues_k_0.rds")
bb_sd_vec <- c(0, 0.4, 0.8)

make_scatter_label_panel <- function(tag, eta) {
  ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0("(", tag, ") "),
      family = style$base_family,
      hjust = 1,
      size = style$geom_text_size_mm(x_axis_title_pt)
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0("eta==", eta),
      parse = TRUE,
      family = style$math_family,
      hjust = 0,
      size = style$geom_text_size_mm(x_axis_title_pt)
    )
}

panel_with_bottom_label <- function(main_plot, label_plot) {
  main_plot / label_plot + plot_layout(heights = c(1, 0.14))
}

scatter_panels <- Map(function(bb_sd, tag) {
  bb_sd_tag <- gsub("\\.", "_", as.character(bb_sd))
  q_name <- paste0("Q_bb_sd_", bb_sd_tag)
  prp_name <- paste0("PRP_bb_sd_", bb_sd_tag)

  if (!(q_name %in% names(pvalue_results)) || !(prp_name %in% names(pvalue_results))) {
    stop("Missing field(s): ", q_name, " or ", prp_name)
  }

  dataa <- data.frame(
    Q = pvalue_results[[q_name]],
    PRP = pvalue_results[[prp_name]]
  )
  if (nrow(dataa) > 1000) {
    dataa <- dataa[sample.int(nrow(dataa), size = 1000, replace = FALSE), , drop = FALSE]
  }

  main_plot <- ggplot(dataa, aes(x = Q, y = PRP)) +
    geom_point(color = "#6186ad", alpha = 0.6, size = 0.00001) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = 1, color = "black", linewidth = 0.25) +
    labs(
      x = "Frequentist p-values",
      y = "Posterior-PRPs"
    ) +
    style$theme +
    theme(
      axis.title = element_text(face = "plain"),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      plot.margin = margin(5, 6, 2, 6)
    )
  panel_with_bottom_label(main_plot, make_scatter_label_panel(tag, bb_sd))
}, bb_sd_vec, c("a", "b", "c"))

scatter_combined <- wrap_plots(scatter_panels, nrow = 1) &
  theme(plot.margin = margin(2, 2, 2, 2))

style$save_pdf(
  scatter_combined,
  "../plot/fixedPRPvsQ_eta_1x3",
  w = 7,
  h = 2.3
)




##############################################################################
######## Plots for checking MCMC convergence
# # Line plot of iteration number vs. bar_beta
##############################################################################
for (i in 1:10) {
  p <- ggplot(
    data.frame(
      iteration = 9000:N,
      bar_beta = bar_betas[i, 9000:N]
    ),
    aes(x = iteration, y = bar_beta)
  ) +
    geom_line() +
    labs(
      title = NULL,
      x = "Iteration",
      y = expression("Sampled" ~ bar(beta))
    )

  print(p)
  Sys.sleep(1) # Pause for visualization
}

hist(bar_betas)

for (i in 1:5) {
  acf(bar_betas[i, (N * r):N])
  Sys.sleep(1) # Pause for visualization
}
