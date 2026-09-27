############################################################################
## Shared header for Fig 1A/C/D/E. source() this before building any panel.
############################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(gridtext)    # italic species names in row labels (gt_render)
library(ggplot2)     # panel E (and B) trajectory plots
library(patchwork)   # combining D + shared legend + E

############################################################################
## 0. Journal page grid (Nature Aging). Width fixed; height is a ceiling
## (warn, don't shrink width, if content goes over).
############################################################################

DOUBLE_COLUMN_WIDTH_MM <- 183
TARGET_WIDTH_MM        <- DOUBLE_COLUMN_WIDTH_MM / 2   # 91.5 mm
MAX_HEIGHT_MM          <- 185

FONT_PT_FLOOR <- 5   # never go smaller than this

############################################################################
## 1. Baseline size constants
############################################################################

base_params <- function() {
  list(
    FONT_PT              = 6,
    ROW_HEIGHT_MM        = 2.6,
    HEATMAP_WIDTH_MM     = 45,   # color body width only
    DEND_WIDTH_MM        = 5,
    LEFT_ANNO_WIDTH_MM   = 2,
    AGE_BAR_HEIGHT_MM    = 2,
    SAMPLE_BAR_HEIGHT_MM = 8,
    LEGEND_WIDTH_MM      = 22,
    LEGEND_GRID_MM       = 2,
    HT_GAP_MM            = 2,
    ROW_GAP_MM           = 0.5
  )
}

## Shared z-score color scale, A/C/D all use the same direction.
Z_LIMIT <- 4
DIVERGING_PALETTE <- c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")   # blue (low) -> red (high)
zscore_col_fun <- colorRamp2(
  c(-Z_LIMIT, -Z_LIMIT / 2, 0, Z_LIMIT / 2, Z_LIMIT),
  DIVERGING_PALETTE
)

## Same palette, different limit -- for panels whose range isn't +-4 (e.g. D).
make_diverging_col_fun <- function(limit) {
  colorRamp2(c(-limit, -limit / 2, 0, limit / 2, limit), DIVERGING_PALETTE)
}

## Colorblind-safe (Okabe-Ito) colors for the 3 feature clusters (D + E).
## C1 is reddish purple, not Okabe-Ito blue -- A/B/C use red for Cytokine
## response and blue for Proteomics, so a blue Cluster 1 in D/E would read
## as "Proteomics" to anyone treating color as consistent across the whole
## figure. Reddish purple keeps the 3-cluster set colorblind-safe while
## staying clear of both that red and that blue.
CLUSTER_COLORS <- c(C1 = "#CC79A7", C2 = "#E69F00", C3 = "#009E73")

## Shared ggplot theme for non-heatmap panels (B, E, ...).
common_ggplot_theme <- function(font_pt) {
  ggplot2::theme_classic(base_size = font_pt) +
    ggplot2::theme(
      axis.text          = ggplot2::element_text(size = font_pt, colour = "black"),
      axis.title         = ggplot2::element_text(size = font_pt),
      axis.line          = ggplot2::element_line(linewidth = 0.3),
      axis.ticks         = ggplot2::element_line(linewidth = 0.3),
      legend.key         = ggplot2::element_blank(),
      legend.background  = ggplot2::element_blank(),
      legend.title       = ggplot2::element_text(size = font_pt, face = "bold"),
      legend.text        = ggplot2::element_text(size = font_pt),
      plot.margin        = ggplot2::margin(1, 1, 1, 1, "mm")
    )
}

## Legend text for every z-score heatmap -- do not reuse as a Heatmap()
## `name=` on more than one heatmap (must be unique per heatmap object).
ZSCORE_LEGEND_NAME <- "Row z-score"

## Age axis / sample bar, shared scale across cohorts.
AGE_BREAKS       <- seq(20, 80, by = 10)
SAMPLE_Y_MAX     <- 80
AGE_COLOR_DOMAIN <- c(18, 85)
age_col_fun <- colorRamp2(AGE_COLOR_DOMAIN, c("#f2e8cf", "#a7c957"))

make_age_labels <- function(age_sorted) {
  ifelse(age_sorted %in% AGE_BREAKS, as.character(age_sorted), "")
}

make_sample_count_anno <- function(sample_counts, p) {
  HeatmapAnnotation(
    `Sample size` = anno_barplot(
      sample_counts,
      gp = gpar(fill = "#a8dadc", col = NA, lwd = 0),
      height = unit(p$SAMPLE_BAR_HEIGHT_MM, "mm"),
      ylim = c(0, SAMPLE_Y_MAX),
      axis_param = list(at = c(0, SAMPLE_Y_MAX / 2, SAMPLE_Y_MAX),
                         gp = gpar(fontsize = p$FONT_PT)),
      border = FALSE
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = p$FONT_PT)
  )
}

## ASCII by default -- Unicode (gamma, middle-dot) caused mbcsToSbcs
## warnings on non-UTF8 devices. Set USE_UNICODE <- TRUE if confirmed safe.
USE_UNICODE <- FALSE
SEP   <- if (USE_UNICODE) " · " else " - "
GAMMA <- if (USE_UNICODE) "γ"   else "gamma"

clean_cytokine_label <- function(x) {
  x <- gsub("C\\.albicanshyphae",   "*C. albicans* hyphae",   x)
  x <- gsub("C\\.albicansconidia",  "*C. albicans* conidia",  x)
  x <- gsub("A\\.fumigatusconidia", "*A. fumigatus* conidia", x)
  x <- gsub("B\\.burgdorferimix",   "*B. burgdorferi* mix",   x)
  x <- gsub("B\\.burgdorferi",      "*B. burgdorferi*",       x)
  x <- gsub("S\\.aureus",           "*S. aureus*",             x)
  x <- gsub("Bacteroides",          "*Bacteroides* spp.",      x)
  x <- sub("^IL6_",  paste0("IL-6", SEP),  x)
  x <- sub("^IL17_", paste0("IL-17", SEP), x)
  x <- sub("^IL22_", paste0("IL-22", SEP), x)
  x <- gsub("_", SEP, x, fixed = TRUE)
  x <- gsub("7days", "7 d",  x)
  x <- gsub("24h",   "24 h", x)
  x <- gsub("48h",   "48 h", x)
  x
}

## Fig 1C: "IL18_Inflammation" -> display "IL18" + category "Inflammation"
split_protein_category <- function(x) {
  category <- sub("^.*_", "", x)
  display  <- sub("_[^_]+$", "", x)
  display  <- gsub("_", "-", display)
  list(display = display, category = category)
}

## Renders gt_render markdown (e.g. "*C. albicans*" -> italics).
render_labels <- function(labels, fontsize) {
  gt_render(labels, gp = gpar(fontsize = fontsize))
}

open_pdf_device <- function(file, width, height) {
  if (USE_UNICODE && capabilities("cairo")) {
    grDevices::cairo_pdf(file, width = width, height = height, onefile = FALSE)
  } else {
    grDevices::pdf(file, width = width, height = height, onefile = FALSE)
  }
}

CAPTION_FONT_PT   <- 5
CAPTION_HEIGHT_MM <- 10

## Selection-criteria caption -- fill in real thresholds before submission.
draw_selection_caption <- function() {
  grid.text(
    "Individuals ordered by age. Features: Spearman age association,\nFDR < 0.05 (fill in actual threshold used), sex-adjusted.",
    x = unit(0.5, "npc"), y = unit(2, "mm"), just = "bottom",
    gp = gpar(fontsize = CAPTION_FONT_PT, col = "grey30")
  )
}

############################################################################
## 2. Fit to the journal grid
############################################################################

## Draws draw_call() on a throwaway device and reads back its true size.
measure_heatmap_size_in <- function(draw_call, max_width_in = 30, max_height_in = 40) {
  grDevices::pdf(NULL, width = max_width_in, height = max_height_in)
  ht_drawn <- draw_call()
  w_in <- grid::convertWidth(ComplexHeatmap:::width(ht_drawn),   "inch", valueOnly = TRUE)
  h_in <- grid::convertHeight(ComplexHeatmap:::height(ht_drawn), "inch", valueOnly = TRUE)
  grDevices::dev.off()
  c(width_in = unname(w_in), height_in = unname(h_in))
}
## For a Legend/Legends object (from Legend() or packLegend()) instead of
## a Heatmap. Can't reuse measure_heatmap_size_in() here -- draw() on a
## Legends object returns NULL, not the object itself, so width()/height()
## must be called on the object directly rather than on draw()'s result.
measure_legend_size_in <- function(lgd, max_width_in = 30, max_height_in = 40) {
  grDevices::pdf(NULL, width = max_width_in, height = max_height_in)
  w_in <- grid::convertWidth(ComplexHeatmap:::width(lgd),   "inch", valueOnly = TRUE)
  h_in <- grid::convertHeight(ComplexHeatmap:::height(lgd), "inch", valueOnly = TRUE)
  grDevices::dev.off()
  c(width_in = unname(w_in), height_in = unname(h_in))
}

## Fields that affect only height -- height_scale multiplies just these.
VERTICAL_ONLY_FIELDS <- c("ROW_HEIGHT_MM", "AGE_BAR_HEIGHT_MM", "SAMPLE_BAR_HEIGHT_MM",
                          "HT_GAP_MM", "ROW_GAP_MM")

## Too wide (long labels) -> shrink everything to fit.
## Too narrow (labels short/hidden) -> widen only the color body.
fit_to_journal_grid <- function(params, build_fn, height_scale = 1, safety = 0.99) {
  p1 <- params
  for (nm in VERTICAL_ONLY_FIELDS) p1[[nm]] <- p1[[nm]] * height_scale

  draw_call_1 <- build_fn(p1)
  size1 <- measure_heatmap_size_in(draw_call_1)
  w1_mm <- size1[["width_in"]] * 25.4

  if (w1_mm >= TARGET_WIDTH_MM) {
    scale <- (TARGET_WIDTH_MM / w1_mm) * safety
    p_final <- p1
    for (nm in names(p_final)) if (is.numeric(p_final[[nm]])) p_final[[nm]] <- p_final[[nm]] * scale
    if (p_final$FONT_PT < FONT_PT_FLOOR) {
      warning(sprintf(
        "Font would drop to %.2f pt to fit %.1f mm width -- clamped to %.0f pt. Width will end up slightly over target. Shorten row labels or hide them to fix.",
        p_final$FONT_PT, TARGET_WIDTH_MM, FONT_PT_FLOOR
      ))
      p_final$FONT_PT <- FONT_PT_FLOOR
    }
    note <- sprintf("row labels are the width bottleneck, rescaled by %.3f", scale)
  } else {
    p_final <- p1
    p_final$HEATMAP_WIDTH_MM <- p_final$HEATMAP_WIDTH_MM + (TARGET_WIDTH_MM - w1_mm)
    note <- "widened color body to fill remaining width; height untouched"
  }

  draw_call_final <- build_fn(p_final)
  size_final <- measure_heatmap_size_in(draw_call_final)
  w_final_mm <- size_final[["width_in"]] * 25.4
  h_final_mm <- size_final[["height_in"]] * 25.4 + CAPTION_HEIGHT_MM

  if (h_final_mm > MAX_HEIGHT_MM) {
    warning(sprintf(
      "Height %.1f mm is %.1f mm over the %.0f mm max. Reduce labelled rows (anno_mark) or lower height_scale -- don't shrink width. Exporting as-is.",
      h_final_mm, h_final_mm - MAX_HEIGHT_MM, MAX_HEIGHT_MM
    ))
  }

  message(sprintf(
    "height_scale %.3f | %s | result %.1f x %.1f mm | font %.2f pt, row height %.2f mm",
    height_scale, note, w_final_mm, h_final_mm, p_final$FONT_PT, p_final$ROW_HEIGHT_MM
  ))

  list(draw_call = draw_call_final, params = p_final,
       width_mm = w_final_mm, height_mm = h_final_mm)
}

## Exports PDF (vector) + PNG (300 dpi), both at exactly TARGET_WIDTH_MM wide.
PNG_DPI <- 300

export_fitted_figure <- function(file_stem, fit, caption_fn = draw_selection_caption) {
  w_in <- TARGET_WIDTH_MM / 25.4
  h_in <- fit$height_mm / 25.4

  pdf_file <- paste0(file_stem, ".pdf")
  open_pdf_device(pdf_file, width = w_in, height = h_in)
  fit$draw_call()
  if (!is.null(caption_fn)) caption_fn()
  grDevices::dev.off()

  png_file <- paste0(file_stem, ".png")
  grDevices::png(png_file, width = w_in * PNG_DPI, height = h_in * PNG_DPI, res = PNG_DPI)
  fit$draw_call()
  if (!is.null(caption_fn)) caption_fn()
  grDevices::dev.off()

  message(sprintf("%s / %s saved at %.1f x %.1f mm.", pdf_file, png_file, TARGET_WIDTH_MM, fit$height_mm))
}
