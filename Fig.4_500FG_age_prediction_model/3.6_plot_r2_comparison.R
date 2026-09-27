# 3.6 Compare age-prediction R² across single omics layers and the multi-omics model (500FG)

library(tidyverse)

# ---- Paths ----
result_dir <- "results/prediction"
fig_dir    <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load results (from 3.1, 3.3 and 3.5) ----
methylation <- readRDS(file.path(result_dir, "methylation_prediction_pvalue_select1e04_results_spearman_elasticnet.rds"))
olink       <- readRDS(file.path(result_dir, "olink_prediction_elastic_results.rds"))
cytokine    <- readRDS(file.path(result_dir, "cytokine_prediction_elastic_results.rds"))
metabolite  <- readRDS(file.path(result_dir, "metabolite_prediction_results.rds"))
microbiome  <- readRDS(file.path(result_dir, "microbiome_prediction_elastic_results.rds"))
cellcount   <- readRDS(file.path(result_dir, "cellcount_prediction_elastic_results.rds"))
multi_omics <- readRDS(file.path(result_dir, "multi_prediction_results_original_methy_separ_elasticnet.rds"))

# ---- Extract R² ----
extract_r2_values <- function(results, data_type) {
  train_r2 <- sapply(results, function(x) x$train_R2)
  val_r2   <- sapply(results, function(x) x$val_R2)
  data.frame(
    Set = rep(c("Training", "Internal test"), each = length(train_r2)),
    R2 = c(train_r2, val_r2),
    DataType = data_type
  )
}

plot_data_r2 <- rbind(
  extract_r2_values(methylation, "Methylation"),
  extract_r2_values(olink, "Proteomics"),
  extract_r2_values(cytokine, "Cytokine"),
  extract_r2_values(metabolite, "Metabolite"),
  extract_r2_values(microbiome, "Microbiome"),
  extract_r2_values(cellcount, "Cell Count"),
  extract_r2_values(multi_omics, "Multi_omics")
)
plot_data_r2$Set      <- factor(plot_data_r2$Set, levels = unique(plot_data_r2$Set))
plot_data_r2$DataType <- factor(plot_data_r2$DataType, levels = unique(plot_data_r2$DataType))

x_labels <- c(
  "Methylation" = "DNA methylation",
  "Proteomics"  = "Proteomics",
  "Cytokine"    = "Cytokine\nresponses",
  "Metabolite"  = "Metabolomics",
  "Microbiome"  = "Microbiomes",
  "Cell Count"  = "Immune\ncell counts",
  "Multi_omics" = "Multi-omics"
)

# ---- Plot (Nature style, 5–6 pt sans-serif) ----
col_training   <- "#7EB6D9"
col_validation <- "#82C09A"

theme_nature <- function() {
  theme_classic(base_family = "sans") +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_blank()
    )
}

g <- ggplot(plot_data_r2, aes(x = DataType, y = R2, fill = Set, color = Set)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75),
               color = "black", alpha = 0.7, linewidth = 0.35) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 0.8, alpha = 0.7) +
  labs(x = NULL, y = expression(R^2)) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = c("Training" = col_training, "Internal test" = col_validation)) +
  scale_color_manual(values = c("Training" = col_training, "Internal test" = col_validation)) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.y = element_text(size = 5, colour = "black"),
    legend.text = element_text(size = 5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    legend.box.spacing = unit(0, "cm"),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# 180 mm x 61 mm (double column)
ggsave(file.path(fig_dir, "3.6_multi-omics_r2_comparison.pdf"), g,
       width = 180 / 25.4, height = 61 / 25.4, device = "pdf", useDingbats = FALSE)
print(g)
