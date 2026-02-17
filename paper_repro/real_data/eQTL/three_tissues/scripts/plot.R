################################################################################
######## Histogram + smoothScatter + forestplots: combined 2x2 figure
################################################################################
rm(list = ls())

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(rstudioapi)

if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

source("../../../../plot_style.R")

# Match the visual style used in batch_effect/scripts/plot_batch.R.
style <- plot_fullwidth_style(
  base_size = 7.5,
  width = 7,
  height = 6.2,
  grid = TRUE
)
x_axis_title_pt <- style$base_size * 1.1

prp <- readRDS("../output/results_PRP.rds")
required_prp_cols <- c("gene", "p_values")
if (!all(required_prp_cols %in% names(prp))) {
  stop(
    "Missing required fields in PRP results: ",
    paste(setdiff(required_prp_cols, names(prp)), collapse = ", ")
  )
}
prp_values <- prp$p_values

eqtlbma <- read.table("../eqtlbma/marginal_probability.txt.gz", header = TRUE)
required_eqtl_cols <- c("gene", "gene.post", "marg.config.1.2.3")
if (!all(required_eqtl_cols %in% names(eqtlbma))) {
  stop(
    "Missing required fields in eqtlbma file: ",
    paste(setdiff(required_eqtl_cols, names(eqtlbma)), collapse = ", ")
  )
}

idx <- match(prp$gene, eqtlbma$gene)
if (anyNA(idx)) {
  stop("Some genes in PRP results are missing in eqtlbma marginal_probability.txt.gz")
}
eqtlbma <- eqtlbma[idx, , drop = FALSE]
if (!identical(prp$gene, eqtlbma$gene)) {
  stop("Gene order mismatch between PRP and eqtlbma after matching")
}

# Tissue-consistency probability from eQtlBma.
eqtl_consistency <- 1 - eqtlbma$gene.post + eqtlbma$marg.config.1.2.3

eqtl_data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
eqtl_data <- na.omit(eqtl_data)

selected_genes <- c(
  "ENSG00000112877.7",
  "ENSG00000123200.16"
)
gene_symbols <- c(
  "CEP72",
  "ZC3H13"
)

make_bottom_panel <- function(axis_title, panel_tag) {
  ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate(
      "text",
      x = 0.5,
      y = 0.62,
      label = axis_title,
      family = style$base_family,
      hjust = 0.5,
      size = style$geom_text_size_mm(x_axis_title_pt)
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.16,
      label = panel_tag,
      family = style$base_family,
      hjust = 0.5,
      size = style$geom_text_size_mm(x_axis_title_pt)
    )
}

panel_with_bottom_label <- function(main_plot, axis_title, panel_tag, label_height = 0.17) {
  main_plot / make_bottom_panel(axis_title, panel_tag) +
    plot_layout(heights = c(1, label_height))
}

hist_bins <- 20
hist_breaks <- seq(0, 1, length.out = hist_bins + 1)
hist_max_count <- max(hist(prp_values, breaks = hist_breaks, plot = FALSE)$counts)
hist_ymax <- max(50, ceiling(hist_max_count * 1.05 / 50) * 50)
hist_ymin <- -0.08 * hist_ymax
mean_prp <- mean(prp_values)

hist_plot <- ggplot(data.frame(p_value = prp_values), aes(x = p_value)) +
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
  labs(
    x = NULL,
    y = "Count"
  ) +
  style$theme +
  theme(
    axis.title = element_text(face = "plain"),
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.margin = margin(5, 6, 2, 6)
  )

my_colramp <- colorRampPalette(brewer.pal(9, "YlGnBu"))

scatter_df <- data.frame(
  posterior_prp = prp_values,
  tissue_consistency = eqtl_consistency
)
if (nrow(scatter_df) > 6000) {
  set.seed(1)
  scatter_points <- scatter_df[sample.int(nrow(scatter_df), 6000, replace = FALSE), , drop = FALSE]
} else {
  scatter_points <- scatter_df
}

smooth_scatter_plot <- ggplot(scatter_df, aes(x = posterior_prp, y = tissue_consistency)) +
  stat_density_2d(
    aes(fill = after_stat(ndensity)),
    geom = "raster",
    contour = FALSE,
    n = 220
  ) +
  geom_point(
    data = scatter_points,
    color = grDevices::adjustcolor("#1f4e79", alpha.f = 0.18),
    shape = 16,
    size = 0.22,
    stroke = 0
  ) +
  scale_fill_gradientn(colors = my_colramp(256), guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0.05)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(xlim = c(0, 0.8183), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = NULL,
    y = "Tissue-consistency probability (eQtlBma)"
  ) +
  style$theme +
  theme(
    axis.title = element_text(face = "plain"),
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.margin = margin(5, 6, 2, 6)
  )

build_forest_plot <- function(gene_id, gene_symbol, clip_range, xticks) {
  selected_gene <- eqtl_data[eqtl_data$gene == gene_id, , drop = FALSE]
  if (nrow(selected_gene) != 1) {
    stop("Cannot find exactly one row for gene: ", gene_id)
  }

  tissues <- c("Artery aorta", "Liver", "Muscle skeletal")
  means <- c(
    selected_gene$bhat_Artery_Aorta,
    selected_gene$bhat_Liver,
    selected_gene$bhat_Muscle_Skeletal
  )
  ses <- c(
    selected_gene$se_Artery_Aorta,
    selected_gene$se_Liver,
    selected_gene$se_Muscle_Skeletal
  )
  lowers <- means - 1.96 * ses
  uppers <- means + 1.96 * ses

  df <- data.frame(
    tissue = factor(tissues, levels = rev(tissues)),
    mean = means,
    lower = lowers,
    upper = uppers
  )

  label_x <- clip_range[2] * 0.35

  ggplot(df, aes(x = mean, y = tissue)) +
    geom_vline(
      xintercept = 0, color = "gray50",
      linetype = "dashed", linewidth = 0.3
    ) +
    geom_segment(
      aes(x = lower, xend = upper, yend = tissue),
      color = "black", linewidth = 0.8
    ) +
    geom_point(shape = 15, size = 2, color = "black") +
    geom_text(
      aes(label = tissue),
      x = label_x, hjust = 0,
      size = style$geom_text_size_mm(style$base_size * 1.1),
      family = style$base_family
    ) +
    scale_x_continuous(breaks = xticks) +
    coord_cartesian(xlim = clip_range) +
    labs(x = NULL, y = NULL) +
    style$theme +
    theme(
      axis.title.x = element_text(
        face = "plain",
        size = x_axis_title_pt
      ),
      axis.text.x = element_text(
        size = style$base_size * 1.0
      ),
      axis.text.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5, 6, 2, 6)
    )
}

clip_ranges <- list(
  c(-1.8, 1.8),
  c(-0.5, 0.5)
)
xticks_list <- list(
  c(-1.8, -1.2, -0.6, 0, 0.6, 1.2, 1.8),
  seq(-0.5, 0.5, by = 0.1)
)

forest_plot_1 <- build_forest_plot(
  selected_genes[1],
  gene_symbols[1],
  clip_ranges[[1]],
  xticks_list[[1]]
)
forest_plot_2 <- build_forest_plot(
  selected_genes[2],
  gene_symbols[2],
  clip_ranges[[2]],
  xticks_list[[2]]
)

panel_a <- panel_with_bottom_label(
  hist_plot,
  "Posterior-PRPs",
  "(a) Histogram of posterior-PRPs"
)
panel_b <- panel_with_bottom_label(
  smooth_scatter_plot,
  "Posterior-PRPs",
  "(b) Smoothed scatter plot"
)
panel_c <- panel_with_bottom_label(
  forest_plot_1,
  paste("Effect size for gene", gene_symbols[1]),
  "(c) Posterior-PRP flagged only"
)
panel_d <- panel_with_bottom_label(
  forest_plot_2,
  paste("Effect size for gene", gene_symbols[2]),
  "(d) Posterior-PRP flagged only"
)

combined_plot <- (panel_a | panel_b) / (panel_c | panel_d) &
  theme(plot.margin = margin(2, 2, 2, 2))

plot_dir <- "../plot"

style$save_pdf(
  combined_plot,
  file.path(plot_dir, "eqtl_prp_2x2_hist_scatter_forest"),
  w = 7,
  h = 6.2
)

# Keep manuscript figure folder in sync when it exists.
bio_figure_dir <- normalizePath("../../../../../Bioinformatics/figure", winslash = "/", mustWork = FALSE)
if (dir.exists(bio_figure_dir)) {
  file.copy(
    file.path(plot_dir, "eqtl_prp_2x2_hist_scatter_forest.pdf"),
    file.path(bio_figure_dir, "eqtl_prp_2x2_hist_scatter_forest.pdf"),
    overwrite = TRUE
  )
}
