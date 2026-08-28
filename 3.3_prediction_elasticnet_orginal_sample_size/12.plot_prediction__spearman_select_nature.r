library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)
library(reshape2)
library(caret)
library(RColorBrewer)

# methylation <- readRDS( "/vol/projects/yzhang/500FG_aging/output/12_prediction/methylation_5000+clock_prediction_results_thrhigh.rds")
# methylation <-readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/methylation_prediction_pvalue_common_select1e04_results_spearman.rds")  
# olink <-readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/olink_prediction_pvalue_common_select_results_spearman.rds")
# cytokine <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/cytokine_prediction_pvalue_common_select_results_spearman.rds")
# metabolite <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/metabolite_prediction_pvalue_common_select_results_spearman.rds")
# microbiome <- readRDS( "/vol/projects/yzhang/500FG_aging/output/12_prediction/microbiome_prediction_pvalue_common_select_results_spearman.rds")
# cellcount <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/cellcounts_prediction_pvalue_common_select_results_spearman.rds")
# Muti_omics <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/spearman/omics_prediction_results_spearman_select.rds")
methylation <-readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/gender/methylation_prediction_pvalue_common_select1e04_results_spearman_gender.rds")  
olink <-readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/olink_prediction_elastic_results.rds")
cytokine <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/cytokine_prediction_elastic_results.rds")
metabolite <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/metabolite_prediction_results.rds")
microbiome <- readRDS( "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/microbiome_prediction_elastic_results.rds")
cellcount <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/cellcount_prediction_elastic_results.rds")
Muti_omics <- readRDS("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/multi_prediction_results_original_methy_separ_elasticnet.rds")
##may cover the originial one

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
plot_data_r2 <- rbind( methylation_r2, olink_r2, cytokine_r2,
                      metabolite_r2, microbiome_r2, 
                      cellcount_r2, Muti_omics_r2)

# Check the combined dataset
head(plot_data_r2)

# Define image resolution and file path
ppi <- 300
fig_width_mm <- 180   # double column width (mm)
fig_height_mm <- 61   # height (mm), adjust if needed
fig_width_in <- fig_width_mm / 25.4
fig_height_in <- fig_height_mm / 25.4

# Save as PDF (vector, preferred for submission) and PNG (300 dpi)
pdf("/vol/projects/yzhang/500FG_aging/output/12_prediction/multi-omics_layers_r2_comparison_orign_sample_nature.pdf",
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