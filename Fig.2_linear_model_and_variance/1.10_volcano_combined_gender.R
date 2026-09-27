# Combined volcano-plot figures for the age~Gender association results
# from all 9 molecular layers (the output of scripts/1.1_..1.10_*_gender.R).
# Run from the repository root. No data is included in this repo -- see
# README for the expected input layout.

library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(grid)
library(ggsci)
library(ggtext)
library(gtable)
library(ggrastr)
library(stringr)
library(tidyverse)

source("utils.R")

combined_dir <- file.path(project_dir, "output", "13_volcano_combined")
platelet_dir <- file.path(project_dir, "output", "09_platelet")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

# age not scaled; Gender-adjusted results from each layer script
olink <- read.csv(file.path(project_dir, "output", "02_Olink", "olink_age_liner_gender_FDR.csv"), row.names = 1)
cytokine <- read.csv(file.path(project_dir, "output", "01_cytokine", "cytokine_age_liner_gender_FDR.csv"), row.names = 1)
metabolite <- read.csv(file.path(project_dir, "output", "03_metabolite", "metabolite_age_linear_gender_FDR.csv"), row.names = 1)
methylation <- readRDS(file.path(project_dir, "output", "05_methylation", "methylation_age_liner_gender_FDR.RDS"))
hormone <- read.csv(file.path(project_dir, "output", "06_hormone", "hormone_age_liner_gender_FDR.csv"), row.names = 1)
cellcounts <- read.csv(file.path(project_dir, "output", "07_cellcounts", "cellcounts_age_rawless_gender_FDR.csv"), row.names = 1)
immunoglobulin <- read.csv(file.path(project_dir, "output", "08_immunoglobulin", "immunoglobulin_age_gender_FDR.csv"), row.names = 1)
microbiome <- read.csv(file.path(project_dir, "output", "10_microbiome", "microbiome_age_gender_FDR.csv"), row.names = 1)
platelet <- read.csv(file.path(project_dir, "output", "09_platelet", "platelet_age_gender_FDR.csv"), row.names = 1)

# Clean up microbiome pathway IDs into readable labels (used for the
# in-plot feature labels)
microbiome <- microbiome %>% rename(feature_old = feature)
microbiome$feature <- microbiome$feature_old %>%
  str_remove("^PWY\\.\\d+\\.\\.?") %>%
  str_remove("^[A-Z0-9]+\\.PWY\\.\\.?") %>%
  str_replace_all("\\.", " ") %>%
  str_squish()

###combined volcano panel, shared x/y axis across layers ###
create_combined_volcano_plot_nature <- function(data_list, titles, output_path) {
  # Title box background color per panel (same order as titles)
  color_palette <- setNames(
  c("#8DD3C7", "#FB8072", "#E9C46A", "#80B1D3", "#BEBADA", "#FCCDE5", "#D98C6B", "#5FAF9D"),
  titles
)

  # Unified colors for all panels: Positive = pink, Negative = blue
  col_pos <- "#FFB3E6"
  col_neg <- "#80B1D3"
  col_ns <- "#B8B8B8"

  # Global x/y limits shared by all panels (padj == 0 would give Inf without pmax)
  limits_df <- bind_rows(data_list) %>%
    mutate(log_padj = -log10(pmax(padj, .Machine$double.xmin)))
  x_range <- range(limits_df$estimate, na.rm = TRUE)
  y_range <- range(limits_df$log_padj, na.rm = TRUE)

  x_span <- diff(x_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  x_range2 <- c(x_range[1] - 0.01 * x_span, x_range[2] + 0.01 * x_span)

  y_hi <- y_range[2] + 0.06 * max(diff(y_range), y_range[2], na.rm = TRUE)
  y_lo <- -0.03 * max(y_range[2], 0.1, na.rm = TRUE)

  y_breaks <- pretty(c(0, y_hi))
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_hi * 1.02]

  plots <- list()
  title_texts <- character(length(data_list))
  title_colors <- character(length(data_list))

  for (i in seq_along(data_list)) {
    data <- data_list[[i]]
    data <- data %>%
      mutate(
        sig = ifelse(padj < 0.05 & estimate > 0, "upregulated",
                     ifelse(padj < 0.05 & estimate < 0, "downregulated", "not significant")),
        log_padj = -log10(padj)
      )

    up_count <- sum(data$sig == "upregulated")
    down_count <- sum(data$sig == "downregulated")
    total_count <- nrow(data)
    up_pct <- round(up_count / total_count * 100, 1)
    down_pct <- round(down_count / total_count * 100, 1)

    top_upregulated <- data %>%
      filter(sig == "upregulated") %>%
      arrange(desc(log_padj)) %>%
      head(3)

    top_downregulated <- data %>%
      filter(sig == "downregulated") %>%
      arrange(desc(log_padj)) %>%
      head(3)

    top_combined <- bind_rows(top_upregulated, top_downregulated)

    category_with_stats <- paste0(
      titles[i], "\n",
      "Neg: ", down_count, " (", down_pct, "%) | ",
      "Pos: ", up_count, " (", up_pct, "%)"
    )
    title_texts[i] <- category_with_stats
    title_colors[i] <- color_palette[titles[i]]

    p <- ggplot(data, aes(x = estimate, y = log_padj)) +
      ggrastr::geom_point_rast(aes(color = sig), size = 1.0, alpha = 0.8, raster.dpi = 300) +
      scale_color_manual(
        values = c(
          "upregulated" = col_pos,
          "downregulated" = col_neg,
          "not significant" = col_ns
        ),
        labels = c(
          "upregulated" = "Positive",
          "downregulated" = "Negative",
          "not significant" = "NS"
        )
      ) +
      geom_text_repel(
        data = top_combined,
        aes(x = estimate, y = log_padj, label = feature),
        size = 1.5,
        color = "black",
        box.padding = 1.1,
        point.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf,
        force = 9,
        segment.size = 0.3,
        segment.color = "grey50"
      ) +
      ggtitle(category_with_stats) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 5, face = "bold", hjust = 0.5, lineheight = 1.1),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5, color = "black"),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.ticks = element_line(size = 0.5, color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, 0, 0, 0),
        legend.spacing.x = unit(4, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = margin(6, 6, 8, 8)
      ) +
      xlab("Estimate") +
      ylab(expression(-log[10](P[adj]))) +
      coord_cartesian(xlim = x_range2, ylim = c(y_lo, y_hi), expand = FALSE) +
      scale_y_continuous(breaks = y_breaks, expand = c(0, 0))
    plots[[i]] <- p
  }

  # Extract single legend (from first plot), then remove legend from all panels
  legend_grob <- cowplot::get_legend(plots[[1]])
  plots_no_leg <- lapply(plots, function(p) p + theme(legend.position = "none"))

  # Convert to grob, then replace title with custom grob (rect + text) so background color shows in PDF/PNG
  plot_grobs <- lapply(plots_no_leg, ggplotGrob)
  n_panels <- length(plot_grobs)
  for (i in seq_len(n_panels)) {
    gt <- plot_grobs[[i]]
    title_idx <- which(gt$layout$name == "title")
    if (length(title_idx) > 0L) {
      title_grob <- gTree(children = gList(
        rectGrob(gp = gpar(fill = title_colors[i], col = "black", lwd = 1)),
        textGrob(
          title_texts[i], x = 0.5, y = 0.5, default.units = "npc",
          gp = gpar(fontsize = 5, fontface = "bold"), just = c("center", "center")
        )
      ))
      gt$grobs[[title_idx]] <- title_grob
      plot_grobs[[i]] <- gt
    }
  }
  combined_grob <- do.call(arrangeGrob, c(plot_grobs, list(ncol = 4, nrow = 2)))

  # Single legend at bottom center; minimal gap above legend
  final_grob <- arrangeGrob(combined_grob, legend_grob, nrow = 2, heights = c(1, 0.055))

  # Square-ish panels: 400mm wide (100mm per panel); height so each volcano ~square + title bar + one legend
  fig_width_mm <- 180
  fig_height_mm <- 200 * 185 / 400

  ggsave(output_path, plot = final_grob, width = fig_width_mm, height = fig_height_mm, units = "mm", device = "pdf")

  print(paste("PDF saved to:", output_path))

  invisible(final_grob)
}

data_list <- list(
  methylation, olink, metabolite, cytokine, cellcounts, microbiome, hormone, immunoglobulin
)
titles <- c(
  "DNA methylation", "Proteomics", "Metabolomics", "Cytokine response",
  "Immune cell counts", "Microbiomes", "Circulating endocrine traits", "Circulating immune traits"
)

create_combined_volcano_plot_nature(
  data_list = data_list,
  titles = titles,
  output_path = file.path(combined_dir, "volcano_plots_gender_nature_small_same_axis.pdf")
)

### Single-panel volcano (same style), for platelet sup fig ###
create_single_volcano <- function(data, title, output_path, title_col = "#A6CEE3") {
  if (!"feature" %in% names(data)) data$feature <- rownames(data)
  data <- data %>% mutate(
    sig = ifelse(padj < 0.05 & estimate > 0, "upregulated",
                 ifelse(padj < 0.05 & estimate < 0, "downregulated", "not significant")),
    log_padj = -log10(padj)
  )
  up_count <- sum(data$sig == "upregulated")
  down_count <- sum(data$sig == "downregulated")
  total_count <- nrow(data)
  title_text <- paste0(
    title, "\nNeg: ", down_count, " (", round(down_count / total_count * 100, 1),
    "%) | Pos: ", up_count, " (", round(up_count / total_count * 100, 1), "%)"
  )
  top_combined <- bind_rows(
    data %>% filter(sig == "upregulated") %>% arrange(desc(log_padj)) %>% head(3),
    data %>% filter(sig == "downregulated") %>% arrange(desc(log_padj)) %>% head(3)
  )
  p <- ggplot(data, aes(x = estimate, y = log_padj)) +
    geom_point(aes(color = sig), size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("upregulated" = "#f1c0e8", "downregulated" = "#90dbf4", "not significant" = "#B8B8B8"),
      labels = c("upregulated" = "Positive", "downregulated" = "Negative", "not significant" = "NS")
    ) +
    geom_text_repel(
      data = top_combined, aes(label = feature), size = 1.4, color = "black",
      box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
      max.overlaps = 20, force = 3, segment.size = 0.3, segment.color = "grey50"
    ) +
    ggtitle(title_text) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 5, face = "bold", hjust = 0.5, lineheight = 1.1),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5, color = "black"),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.ticks = element_line(size = 0.5, color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0),
      legend.spacing.x = unit(4, "pt"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(4, 6, 6, 4)
    ) +
    xlab("Estimate") + ylab(expression(-log[10](P[adj])))
  gt <- ggplotGrob(p + theme(legend.position = "none"))
  title_idx <- which(gt$layout$name == "title")
  if (length(title_idx) > 0L) {
    gt$grobs[[title_idx]] <- gTree(children = gList(
      rectGrob(gp = gpar(fill = title_col, col = "black", lwd = 1)),
      textGrob(title_text, x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = 5, fontface = "bold"), just = c("center", "center"))
    ))
  }
  leg <- get_legend(p)
  out <- arrangeGrob(gt, leg, nrow = 2, heights = c(1, 0.055))
  w_mm <- 60
  h_mm <- 60 # same aspect as combined
  ggsave(output_path, plot = out, width = w_mm, height = h_mm, units = "mm", device = "pdf")
  ggsave(sub("\\.pdf$", ".png", output_path, ignore.case = TRUE), plot = out, width = w_mm, height = h_mm, units = "mm", dpi = 300, device = "png")
  invisible(out)
}

create_single_volcano(
  platelet, "Platelet",
  file.path(platelet_dir, "platelet_volcano_gender.pdf")
)
