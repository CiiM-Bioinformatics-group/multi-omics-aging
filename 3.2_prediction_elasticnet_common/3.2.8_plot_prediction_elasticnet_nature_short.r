library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)
library(reshape2)
library(caret)
library(RColorBrewer)


methylation <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/methylation_prediction_pvalue_common_results_spearman_elasticnet.rds")
olink <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/olink_prediction_pvalue_common_select_results_spearman_quoteID_elasticnet.rds")
cytokine <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/cytokine_prediction_pvalue_common_results_spearman_quoteID_elasticnet.rds")
metabolite <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/metabolite_prediction_pvalue_common_select_results_spearman_quoteID_elasticnet.rds")
microbiome <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/microbiome_prediction_pvalue_common_select_results_spearman_elasticnet.rds")
cellcount <- readRDS( "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/cellcounts_prediction_pvalue_common_select_results_spearman_elasticnet.rds")
Muti_omics <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/multi_prediction_results_original_methy_separ_elasticnet.rds")


# Extract R² values from each dataset
extract_r2_values <- function(results, data_type) {
  train_r2_values <- sapply(results, function(x) x$train_R2)
  val_r2_values <- sapply(results, function(x) x$val_R2)

# Create a combined dataframe
  data.frame(
    Set = rep(c("Training", "Validation"), each = length(train_r2_values)),
    R2 = c(train_r2_values, val_r2_values),
    DataType = data_type
  )
}

# Create R² dataframes for each dataset
#genotype_r2 <- extract_r2_values(genotype, "Genotype")
methylation_r2 <- extract_r2_values(methylation, "Methylation")
olink_r2 <- extract_r2_values(olink, "Proteomics")
cytokine_r2 <- extract_r2_values(cytokine, "Cytokine")
metabolite_r2 <- extract_r2_values(metabolite, "Metabolite")
microbiome_r2 <- extract_r2_values(microbiome, "Microbiome")
cellcount_r2 <- extract_r2_values(cellcount, "Cell Count")
Muti_omics_r2 <- extract_r2_values(Muti_omics, "Muti_omics")

# Combine all dataframes into one
plot_data_r2 <- rbind(methylation_r2, olink_r2, cytokine_r2,
                      metabolite_r2, microbiome_r2, cellcount_r2, Muti_omics_r2)

# Check the combined dataset
head(plot_data_r2)

 # Define image resolution and file path
ppi <- 300
fig_width_mm <- 120   # double column width (mm)
fig_height_mm <- 61   # height (mm), adjust if needed
fig_width_in <- fig_width_mm / 25.4
fig_height_in <- fig_height_mm / 25.4

# Save as PDF (vector, preferred for submission) and PNG (300 dpi)
pdf("/vol/projects/yzhang/500FG_aging/output/12_prediction/multi-omics_layers_r2_comparison_clock_elasticnet_nature.pdf",
    width = fig_width_in, height = fig_height_in, useDingbats = FALSE)
# Alternatively PNG at 300 dpi:
# png("/vol/projects/yzhang/500FG_aging/output/12_prediction/multi-omics_layers_r2_comparison_clock_nogeno.png",
#     width = fig_width_in * ppi, height = fig_height_in * ppi, res = ppi)

# Nature Aging colors (sans-serif friendly, RGB)
col_training  <- "#7EB6D9"   # blue
col_validation <- "#82C09A"  # green

# Ensure the 'Set' order follows the original data order
plot_data_r2$Set <- factor(plot_data_r2$Set, levels = unique(plot_data_r2$Set))
plot_data_r2$DataType <- factor(plot_data_r2$DataType, levels = unique(plot_data_r2$DataType))

x_labels <- c(
  "Methylation"   = "DNA methylation",
  "Proteomics"    = "Proteomics",
  "Cytokine"      = "Cytokine\nresponses",
  "Metabolite"    = "Metabolomics",
  "Microbiome"    = "Microbiomes",
  "Cell Count"    = "Immune\ncell counts",
  "Muti_omics"    = "Multi-omics"
)
                          
# Create the plot (Nature: 5-7 pt font, sans-serif)
g <- ggplot(plot_data_r2, aes(x = DataType, y = R2, fill = Set, color = Set)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75),
               color = "black", alpha = 0.7, linewidth = 0.35) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 0.8, alpha = 0.7) +
  labs(x = "Omics data type",
       y = expression(R^2)) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = c("Training" = col_training, "Validation" = col_validation)) +
  scale_color_manual(values = c("Training" = col_training, "Validation" = col_validation)) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title = element_text(size = 5, colour = "black"),
    legend.text = element_text(size = 5, colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
  )

# Print and save the plot
print(g)
dev.off()
print(g)          


# ============================
# Additional analysis/plot:
# Paired t-test for Validation R2
# between Methylation and Muti_omics
# ============================

# Pair by iteration only
m_df <- tibble(
  iteration = sapply(methylation, function(x) x$iteration),
  val_r2_m = sapply(methylation, function(x) x$val_R2)
)

o_df <- tibble(
  iteration = sapply(Muti_omics, function(x) x$iteration),
  val_r2_o = sapply(Muti_omics, function(x) x$val_R2)
)

# Keep only iterations present in both datasets
pair_df <- inner_join(m_df, o_df, by = "iteration")
pair_df <- pair_df[complete.cases(pair_df$val_r2_m, pair_df$val_r2_o), ]

if (nrow(pair_df) < 2) {
  stop("Not enough paired validation values after matching by iteration (need at least 2 pairs).")
}

methylation_val_common <- pair_df$val_r2_m
multiomics_val_common <- pair_df$val_r2_o

# Paired t-test (Validation only)
paired_t_res <- t.test(methylation_val_common, multiomics_val_common, paired = TRUE)
print(paired_t_res)
message("Paired t-test pairing mode: iteration-based; N pairs = ", nrow(pair_df))

# Replot with original settings and add significance only if p < 0.05
g_validation_sig <- ggplot(plot_data_r2, aes(x = DataType, y = R2, fill = Set, color = Set)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75),
               color = "black", alpha = 0.7, linewidth = 0.35) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 0.8, alpha = 0.7) +
  labs(x = "Omics data type",
       y = expression(R^2)) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = c("Training" = col_training, "Validation" = col_validation)) +
  scale_color_manual(values = c("Training" = col_training, "Validation" = col_validation)) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title = element_text(size = 5, colour = "black"),
    legend.text = element_text(size = 5, colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    legend.box.spacing = unit(0, "cm"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "pt")
  ) +
  # Y-axis ticks/labels stop at 1; bracket/stars drawn above panel (clip off)
  coord_cartesian(ylim = c(NA, 1), clip = "off")

# Significance stars (state thresholds in figure legend). Prism-style fourth star:
# * p <= 0.05, ** p <= 0.01, *** p <= 0.001, **** p <= 0.0001
pval <- paired_t_res$p.value
sig_label <- if (pval <= 0.0001) {
  "****"
} else if (pval <= 0.001) {
  "***"
} else if (pval <= 0.01) {
  "**"
} else if (pval <= 0.05) {
  "*"
} else {
  ""
}

if (nzchar(sig_label)) {
  # Position for the Validation boxes (2nd level in Set with dodge width 0.75)
  dodge_offset_validation <- 0.75 / 4
  x1 <- which(levels(plot_data_r2$DataType) == "Methylation") + dodge_offset_validation
  x2 <- which(levels(plot_data_r2$DataType) == "Muti_omics") + dodge_offset_validation

  # Bracket and stars above R^2 = 1 (line higher than before; star above line)
  y_line <- 1.09
  y_star <- 1.15

  g_validation_sig <- g_validation_sig +
    theme(plot.margin = margin(22, 4, 4, 4, "pt")) +
    annotate("segment", x = x1, xend = x2, y = y_line, yend = y_line,
             linewidth = 0.35, colour = "black") +
    annotate("text", x = (x1 + x2) / 2, y = y_star, label = sig_label,
             size = 3, colour = "black")
}

# Save the new plot
pdf("/vol/projects/yzhang/500FG_aging/output/12_prediction/multi-omics_layers_r2_comparison_clock_elasticnet_nature_validation_paired_ttest.pdf",
    width = fig_width_in, height = fig_height_in, useDingbats = FALSE)
print(g_validation_sig)
dev.off()
print(g_validation_sig)