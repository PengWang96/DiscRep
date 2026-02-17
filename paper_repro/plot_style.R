# Common plotting style for journal figures.
# Defaults: full-width (double-column) and single-column presets, PDF output via cairo + showtext.

# ==============================================================================
# Global Font Path Configuration
# ==============================================================================

#' @title Global Font Path Configuration
#' @description
#' Modify these paths if your fonts are installed in different locations.

FONT_PATHS <- list(
  # English fonts (Arial)
  english = c(
    "C:/Windows/Fonts/arial.ttf"
  ),
  # Chinese fonts (SimSun)
  chinese = c(
    "C:/Windows/Fonts/simsun.ttc",
    "C:/Windows/Fonts/simsun.ttf"
  ),

  # Math fonts (Latin Modern Math)
  math = c(
    "C:/CTEX/MiKTeX/fonts/opentype/public/lm-math/latinmodern-math.otf",
    "C:/texlive/2023/texmf-dist/fonts/opentype/public/lm-math/latinmodern-math.otf",
    "C:/texlive/2024/texmf-dist/fonts/opentype/public/lm-math/latinmodern-math.otf"
  )
)

# Default font family names
FONT_FAMILIES <- list(
  english = "Arial",
  chinese = "SimSun",
  math = "Latin Modern Math"
)

# ==============================================================================
# Legend Position Presets
# ==============================================================================

#' @title Legend Position Presets
#' @description
#' Commonly used legend positions for quick access.
#' Use with ggplot2::theme(legend.position = LEGEND_POS$xxx)

LEGEND_POS <- list(
  # Standard positions
  top = "top",
  bottom = "bottom",
  left = "left",
  right = "right",
  none = "none",

  # Inside plot area positions (x, y coordinates from 0-1)
  top_left = c(0.02, 0.98),
  top_right = c(0.98, 0.98),
  bottom_left = c(0.02, 0.02),
  bottom_right = c(0.98, 0.02),
  center = c(0.5, 0.5),

  # Inside with margin offsets
  inside_top = c(0.5, 0.95),
  inside_bottom = c(0.5, 0.05)
)

#' @title Get Legend Theme for Inside Positions
#' @description
#' Returns a ggplot2 theme object with legend positioned inside the plot.
#'
#' @param position A character string or numeric vector of length 2.
#'   Use LEGEND_POS presets or custom c(x, y) coordinates.
#' @param direction "horizontal" or "vertical" legend layout.
#' @param bg_alpha Background alpha for legend box (0-1).
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' p + legend_inside(LEGEND_POS$top_right)
#' p + legend_inside(c(0.8, 0.2), direction = "horizontal")
legend_inside <- function(position = LEGEND_POS$top_right,
                          direction = "vertical",
                          bg_alpha = 0.8) {
  if (is.character(position) && length(position) == 1) {
    return(ggplot2::theme(
      legend.position = position,
      legend.direction = direction
    ))
  }

  ggplot2::theme(
    legend.position = position,
    legend.justification = position,
    legend.direction = direction,
    legend.background = ggplot2::element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = bg_alpha),
      color = NA
    ),
    legend.key = ggplot2::element_rect(fill = "transparent")
  )
}

# ==============================================================================
# Subplot Label Automation
# ==============================================================================

#' @title Generate Subplot Labels
#' @description
#' Generate a sequence of subplot labels in various formats.
#'
#' @param n Number of labels to generate.
#' @param style Label style: "paren_lower" for (a), (b), (c)...;
#'   "paren_upper" for (A), (B), (C)...;
#'   "lower" for a, b, c...;
#'   "upper" for A, B, C...;
#'   "number" for 1, 2, 3...;
#'   "paren_number" for (1), (2), (3)...
#'
#' @return A character vector of labels.
#'
#' @examples
#' subplot_labels(3)
#' # Returns: c("(a)", "(b)", "(c)")
subplot_labels <- function(n, style = "paren_lower") {
  if (n < 1) {
    return(character(0))
  }
  if (n > 26 && style %in% c("paren_lower", "paren_upper", "lower", "upper")) {
    warning("More than 26 labels requested; cycling through alphabet.")
  }

  labels <- switch(style,
    "paren_lower" = paste0("(", letters[((seq_len(n) - 1) %% 26) + 1], ")"),
    "paren_upper" = paste0("(", LETTERS[((seq_len(n) - 1) %% 26) + 1], ")"),
    "lower" = letters[((seq_len(n) - 1) %% 26) + 1],
    "upper" = LETTERS[((seq_len(n) - 1) %% 26) + 1],
    "number" = as.character(seq_len(n)),
    "paren_number" = paste0("(", seq_len(n), ")"),
    stop("Unknown style: ", style)
  )

  return(labels)
}

#' @title Add Subplot Labels to a List of Plots
#' @description
#' Adds tags to a list of ggplot objects for subplot labeling.
#'
#' @param plots A list of ggplot objects.
#' @param style Label style (see subplot_labels).
#' @param tag_size Font size for tags.
#' @param tag_face Font face for tags ("plain", "bold", "italic", "bold.italic").
#'
#' @return A list of ggplot objects with tags added.
#'
#' @examples
#' plots <- list(p1, p2, p3)
#' labeled_plots <- add_subplot_labels(plots)
add_subplot_labels <- function(plots,
                               style = "paren_lower",
                               tag_size = 14,
                               tag_face = "bold") {
  if (!is.list(plots)) {
    stop("plots must be a list of ggplot objects")
  }

  labels <- subplot_labels(length(plots), style = style)

  mapply(function(p, label) {
    p + ggplot2::labs(tag = label) +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(size = tag_size, face = tag_face)
      )
  }, plots, labels, SIMPLIFY = FALSE)
}

# ==============================================================================
# Font Loading Helper
# ==============================================================================

#' @title Find and Load Font
#' @description
#' Internal function to find and load a font from candidate paths.
#'
#' @param candidates Character vector of candidate font paths.
#' @param family_name Name to register the font under.
#'
#' @return The path to the found font, or NA if not found.
find_and_load_font <- function(candidates, family_name) {
  font_file <- candidates[file.exists(candidates)][1]

  if (!is.na(font_file) && !(family_name %in% sysfonts::font_families())) {
    sysfonts::font_add(family = family_name, regular = font_file)
  }

  return(font_file)
}

# ==============================================================================
# Main Style Functions
# ==============================================================================

#' @title Full-Width (Double-Column) Plot Style
#' @description
#' Creates a plot style configuration suitable for full-width figures in
#' double-column journal layouts. Default width is 7 inches.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family for text (default: Arial).
#' @param math_family Font family for mathematical expressions (default: Latin Modern Math).
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param axis_title_rel Relative size multiplier for axis titles.
#' @param axis_text_rel Relative size multiplier for axis tick labels.
#' @param legend_title_rel Relative size multiplier for legend title.
#' @param legend_text_rel Relative size multiplier for legend text.
#' @param tag_size Font size for subplot tags (a, b, c...).
#' @param tag_face Font face for subplot tags.
#' @param tag_style Style for subplot labels (see subplot_labels).
#' @param grid Whether to show grid lines.
#' @param axis_line_width Width of axis lines.
#'
#' @return A list containing theme, save functions, and font information.
#'
#' @examples
#' style <- plot_fullwidth_style()
#' p <- ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   style$theme
#' style$save_pdf(p, "my_figure")
plot_fullwidth_style <- function(
  base_size = 11,
  base_family = FONT_FAMILIES$english,
  math_family = FONT_FAMILIES$math,
  width = 7,
  height = 4,
  axis_title_rel = 1.1,
  axis_text_rel = 0.9,
  legend_title_rel = 1.0,
  legend_text_rel = 0.9,
  tag_size = 14,
  tag_face = "bold",
  tag_style = "paren_lower",
  grid = FALSE,
  axis_line_width = 0.5
) {
  plot_style_common(
    base_size = base_size,
    base_family = base_family,
    math_family = math_family,
    width = width,
    height = height,
    axis_title_rel = axis_title_rel,
    axis_text_rel = axis_text_rel,
    legend_title_rel = legend_title_rel,
    legend_text_rel = legend_text_rel,
    tag_size = tag_size,
    tag_face = tag_face,
    tag_style = tag_style,
    grid = grid,
    axis_line_width = axis_line_width
  )
}

#' @title Single-Column Plot Style
#' @description
#' Creates a plot style configuration suitable for single-column figures.
#' Default width is 3.4 inches (typical journal single column).
#'
#' @inheritParams plot_fullwidth_style
#'
#' @return A list containing theme, save functions, and font information.
#'
#' @examples
#' style <- plot_singlecolumn_style()
#' p <- ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   style$theme
#' style$save_pdf(p, "my_figure")
plot_singlecolumn_style <- function(
  base_size = 10,
  base_family = FONT_FAMILIES$english,
  math_family = FONT_FAMILIES$math,
  width = 3.4,
  height = 3,
  axis_title_rel = 1.1,
  axis_text_rel = 0.9,
  legend_title_rel = 1.0,
  legend_text_rel = 0.9,
  tag_size = 12,
  tag_face = "bold",
  tag_style = "paren_lower",
  grid = FALSE,
  axis_line_width = 0.5
) {
  plot_style_common(
    base_size = base_size,
    base_family = base_family,
    math_family = math_family,
    width = width,
    height = height,
    axis_title_rel = axis_title_rel,
    axis_text_rel = axis_text_rel,
    legend_title_rel = legend_title_rel,
    legend_text_rel = legend_text_rel,
    tag_size = tag_size,
    tag_face = tag_face,
    tag_style = tag_style,
    grid = grid,
    axis_line_width = axis_line_width
  )
}

#' @title Common Plot Style Implementation
#' @description
#' Internal function that implements the common plot style logic.
#'
#' @inheritParams plot_fullwidth_style
#'
#' @return A list containing:
#'   \item{theme}{A ggplot2 theme object}
#'   \item{base_family}{The base font family used}
#'   \item{math_family}{The math font family used}
#'   \item{base_size}{The base font size}
#'   \item{tag_style}{The subplot label style}
#'   \item{geom_text_size_mm}{Function to convert pt to mm for geom_text}
#'   \item{save_pdf}{Function to save plots as PDF}
#'   \item{labels}{Function to generate subplot labels}
#'   \item{add_labels}{Function to add labels to a list of plots}
#'   \item{legend_inside}{Function to position legend inside plot}
plot_style_common <- function(
  base_size,
  base_family,
  math_family,
  width,
  height,
  axis_title_rel,
  axis_text_rel,
  legend_title_rel,
  legend_text_rel,
  tag_size,
  tag_face,
  tag_style,
  grid,
  axis_line_width
) {
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("sysfonts", quietly = TRUE) ||
    !requireNamespace("showtext", quietly = TRUE)) {
    stop("Required packages: ggplot2, sysfonts, showtext")
  }

  # Load base font by family selection (English default, Chinese optional)
  base_font_candidates <- if (identical(base_family, FONT_FAMILIES$chinese)) {
    FONT_PATHS$chinese
  } else {
    FONT_PATHS$english
  }
  font_file <- find_and_load_font(base_font_candidates, base_family)
  if (is.na(font_file)) {
    warning("Font (", base_family, ") not found. Using system default.")
    base_family <- "sans"
  }

  # Load Math font (Latin Modern Math)
  math_file <- find_and_load_font(FONT_PATHS$math, math_family)
  if (is.na(math_file)) {
    math_family <- base_family
  }

  # Enable showtext
  showtext::showtext_auto(enable = TRUE)

  # Build theme
  theme <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = base_size * axis_title_rel,
        face = "bold",
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = base_size * axis_text_rel,
        color = "black"
      ),
      legend.title = ggplot2::element_text(
        size = base_size * legend_title_rel,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(size = base_size * legend_text_rel),
      panel.grid.major = if (grid) {
        ggplot2::element_line(color = "grey85")
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = if (grid) {
        ggplot2::element_line(color = "grey90")
      } else {
        ggplot2::element_blank()
      },
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = axis_line_width
      ),
      plot.tag = ggplot2::element_text(
        size = tag_size,
        face = tag_face,
        family = base_family
      ),
      plot.tag.position = c(0, 1)
    )

  # Helper function for text size conversion
  geom_text_size_mm <- function(pt = base_size) {
    pt / 2.8
  }

  # PDF save function
  save_pdf <- function(plot_obj, filename, w = width, h = height) {
    if (tools::file_ext(filename) == "") {
      filename <- paste0(filename, ".pdf")
    }
    grDevices::cairo_pdf(filename, width = w, height = h)
    showtext::showtext_begin()
    print(plot_obj)
    showtext::showtext_end()
    grDevices::dev.off()
    message("Saved: ", filename)
  }

  # Return style object
  list(
    theme = theme,
    base_family = base_family,
    math_family = math_family,
    base_size = base_size,
    tag_style = tag_style,
    geom_text_size_mm = geom_text_size_mm,
    save_pdf = save_pdf,

    # Subplot label helpers
    labels = function(n) subplot_labels(n, style = tag_style),
    add_labels = function(plots) {
      add_subplot_labels(plots, style = tag_style, tag_size = tag_size, tag_face = tag_face)
    },

    # Legend positioning
    legend_inside = legend_inside,
    LEGEND_POS = LEGEND_POS
  )
}

