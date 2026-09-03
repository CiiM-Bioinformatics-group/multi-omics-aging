library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(grid)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

# load("/vol/projects/yzhang/500FG_aging/output/18_multi_omics_heatmap/multi_omics_heatmap.Rdata")

methylation <- readRDS("/vol/projects/yzhang/500FG_aging/output/05_methylation/methylation_age_liner_gender_FDR.RDS")
methylation <- methylation %>% rename(feature = Mvalue)
head(methylation)

table(methylation$sig)

methylation_filtered <- methylation[methylation$padj < 0.0001 & 
                                     abs(methylation$estimate) > 0.02, ]
dim(methylation_filtered)

olink <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv",row.names=1)  #age not scale
cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv",row.names=1) #age not scale
metabolite <- read.csv("/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR.csv",row.names=1) #age not scale
hormone <- read.csv("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR.csv",row.names=1) ##age not scale
cellcounts <- read.csv("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv",row.names=1)  ##age not scale
immunoglobulin <- read.csv("/vol/projects/yzhang/500FG_aging/output/08_immunoglobulin/immunoglobulin_age_gender_FDR.csv",row.names=1) ##age not scale
microbiome <- read.csv("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR.csv",row.names=1) ##age not scale

head(olink)
head(cytokine)
head(metabolite)
head(cellcounts)
head(microbiome)
head(hormone)
head(immunoglobulin)

Mvalue_raw <- readRDS("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.rds")
head(Mvalue_raw)
dim(Mvalue_raw)

Mvalue_subset <- Mvalue_raw[, colnames(Mvalue_raw) %in% methylation_filtered$feature]
head(Mvalue_subset)
dim(Mvalue_subset)

Mvalue_raw <- NULL
rm(Mvalue_raw)

olink_raw  <- readRDS("/vol/projects/CIIM/cohorts_old/500FG/olinkData/NPX_FG500.RDS")
head(olink_raw)
dim(olink_raw)
cytokine_raw <- read.csv("/vol/projects/yzhang/500FG_aging/input/cytokine/filtered_cytokine_filled_trans.csv",row.names=1)
head(cytokine_raw)
dim(cytokine_raw)
metabolite_raw <- readRDS("/vol/projects/yzhang/500FG_aging/input/metabolite/metabolite_trans.rds")
head(metabolite_raw)
dim(metabolite_raw)
hormone_raw <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_log2_hormone_levels.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
head(hormone_raw)
dim(hormone_raw)
cellcounts_raw <- read.csv("/vol/projects/yzhang/500FG_aging/input/cellCounts/cellcounts_name_replace.csv",row.names=1,check.names = FALSE)
head(cellcounts_raw)
dim(cellcounts_raw)
microbiome_raw <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_microbiome_pathways.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
head(microbiome_raw)
dim(microbiome_raw)
immunoglobulin_raw <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_log2_immunoglobulin_levels.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
head(immunoglobulin_raw)
dim(immunoglobulin_raw)

##nature style
# ============================================
# 1. Filter significant features from each layer
# ============================================
filter_sig <- function(df, layer, is_filtered = FALSE) {
  if (!is_filtered) df <- df[df$sig == "sig", ]
  data.frame(feature = df$feature, estimate = df$estimate, layer = layer)
}

feature_info <- rbind(
  filter_sig(methylation_filtered, "DNA methylation", TRUE),
  filter_sig(olink, "Proteomics"),
  filter_sig(cytokine, "Cytokine responses"),
  filter_sig(metabolite, "Metabolomics"),
  filter_sig(cellcounts, "Immune cell counts"),
  filter_sig(microbiome, "Microbiomes"),
  filter_sig(hormone, "Circulating endocrine traits"),
  filter_sig(immunoglobulin, "Circulating immune traits")
)

# ============================================
# 2. Select sig features from each raw data
# ============================================
select_features <- function(raw_df, sig_features) {
  cols <- colnames(raw_df)[colnames(raw_df) %in% sig_features]
  raw_df[, cols, drop = FALSE]
}

# Mvalue_sub <- select_features(Mvalue_raw, methylation_filtered$feature)
Mvalue_sub <- Mvalue_subset 
olink_sub <- select_features(olink_raw, feature_info$feature[feature_info$layer == "Proteomics"])
cytokine_sub <- select_features(cytokine_raw, feature_info$feature[feature_info$layer == "Cytokine responses"])
metabolite_sub <- select_features(metabolite_raw, feature_info$feature[feature_info$layer == "Metabolomics"])
cellcounts_sub <- select_features(cellcounts_raw, feature_info$feature[feature_info$layer == "Immune cell counts"])
microbiome_sub <- select_features(microbiome_raw, feature_info$feature[feature_info$layer == "Microbiomes"])
hormone_sub <- select_features(hormone_raw, feature_info$feature[feature_info$layer == "Circulating endocrine traits"])
immunoglobulin_sub <- select_features(immunoglobulin_raw, feature_info$feature[feature_info$layer == "Circulating immune traits"])

# ============================================
# 3. Find common samples across all datasets
# ============================================
common_samples <- Reduce(intersect, list(
  rownames(Mvalue_sub),
  rownames(olink_sub),
  rownames(cytokine_sub),
  rownames(metabolite_sub),
  rownames(cellcounts_sub),
  rownames(microbiome_sub),
  rownames(hormone_sub),
  rownames(immunoglobulin_sub)
))

cat("Number of common samples:", length(common_samples), "\n")

# ============================================
# 4. Subset to common samples and merge
# ============================================
raw_data <- cbind(
  Mvalue_sub[common_samples, ],
  olink_sub[common_samples, ],
  cytokine_sub[common_samples, ],
  metabolite_sub[common_samples, ],
  cellcounts_sub[common_samples, ],
  microbiome_sub[common_samples, ],
  hormone_sub[common_samples, ],
  immunoglobulin_sub[common_samples, ]
)

# Update feature_info to only include features that exist in raw_data
feature_info <- feature_info[feature_info$feature %in% colnames(raw_data), ]
# Fix layer order to match cbind order (so heatmap blocks match merge order)
layer_order <- c(
  "DNA methylation", "Proteomics", "Cytokine responses", "Metabolomics",
  "Immune cell counts", "Microbiomes", "Circulating endocrine traits", "Circulating immune traits"
)
feature_info$layer <- factor(feature_info$layer, levels = layer_order)
# ============================================
# 5. Calculate feature-feature correlation matrix
# ============================================
cor_matrix <- cor(raw_data, use = "pairwise.complete.obs", method = "spearman")

# ============================================
# 6. Prepare annotation information
# ===========================================

# Create layer annotation
layer_annotation <- feature_info$layer
names(layer_annotation) <- feature_info$feature

# Age correlation annotation (based on estimate sign)
age_cor <- ifelse(feature_info$estimate > 0, "Positive", "Negative")
names(age_cor) <- feature_info$feature

# Order by layer (clustering will be done within each layer)
feature_order <- feature_info$feature[order(feature_info$layer)]
cor_matrix_ordered <- cor_matrix[feature_order, feature_order]

# Get p-value for each pair from correlation (t-approximation for Spearman)
n <- nrow(raw_data)
t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
t_stat[abs(cor_matrix) >= 1 | is.na(cor_matrix)] <- 0
p_matrix <- 2 * pt(-abs(t_stat), n - 2)
p_matrix[is.na(cor_matrix)] <- NA

# FDR (BH) correction over unique pairs (upper triangle only)
p_adj_matrix <- matrix(NA, nrow(p_matrix), ncol(p_matrix))
p_vec <- p_matrix[upper.tri(p_matrix)]
p_adj_vec <- p.adjust(p_vec, method = "BH")
p_adj_matrix[upper.tri(p_adj_matrix)] <- p_adj_vec
p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
diag(p_adj_matrix) <- NA
p_adj_matrix[is.na(p_matrix)] <- NA

# Heatmap value: direction from correlation, strength from FDR -> signed -log10(FDR)
signed_logp <- sign(cor_matrix) * (-log10(p_adj_matrix))
# Avoid Inf and cap extreme -log10(p) for better color scale (e.g. p < 1e-10)
signed_logp[is.infinite(signed_logp)] <- NA
cap <- 10  # -log10(1e-10) = 10; adjust if needed
signed_logp[signed_logp > cap] <- cap
signed_logp[signed_logp < -cap] <- -cap

# Matrix for heatmap: signed -log10(p), same order
signed_logp_ordered <- signed_logp[feature_order, feature_order]

# Color scale for -log10(p) (symmetric around 0)
max_val <- max(abs(signed_logp_ordered), na.rm = TRUE)

layer_colors <- c(
  "DNA methylation" = "#B3DE69",
  "Proteomics" = "#FB8072", 
  "Cytokine responses" = "#80B1D3",
  "Metabolomics" = "#FFFFB3",
  "Immune cell counts" = "#BEBADA",
  "Microbiomes" = "#FCCDE5",
  "Circulating endocrine traits" = "#FDB462",
  "Circulating immune traits" = "#8DD3C7"
)
   
#age_colors <- c("Positive" = "#c51b7d", "Negative" = "#4575b4")
age_colors <- c("Positive" = "#FB8072", "Negative" = "#80B1D3")
# ============================================
# 8. Create heatmap
# ============================================
gp_small <- gpar(fontsize = 5)
grid_small <- unit(2, "mm") 
anno_legend_param <- list(
  `Age correlation` = list(
    title_gp = gp_small, labels_gp = gp_small,
    grid_height = grid_small, grid_width = grid_small, nrow = 2
  ),
  Layer = list(
    title_gp = gp_small, labels_gp = gp_small,
    grid_height = grid_small, grid_width = grid_small, nrow = 4
  )
)

row_ha <- rowAnnotation(
  `Age correlation` = age_cor[feature_order],
  Layer = layer_annotation[feature_order],
  col = list(`Age correlation` = age_colors, Layer = layer_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gp_small,
  simple_anno_size = unit(2, "mm")
)

col_ha <- HeatmapAnnotation(
  `Age correlation` = age_cor[feature_order],
  Layer = layer_annotation[feature_order],
  col = list(`Age correlation` = age_colors, Layer = layer_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_name_gp = gp_small,
  annotation_legend_param = anno_legend_param,
  simple_anno_size = unit(2, "mm")
)

col_fun <- colorRamp2(c(-1, 0, 1), c("#669bbc", "white", "#c1121f")) #4575b4 #c51b7d or #80B1D3 #FB8072 ""

ht <- Heatmap(
  signed_logp_ordered,
  name = "-log10 FDR",
  col = col_fun,
  top_annotation = col_ha,
  left_annotation = row_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_split = layer_annotation[feature_order],
  column_split = layer_annotation[feature_order],
  row_title = NULL,
  column_title = NULL,
  border = FALSE,
    row_gap = unit(0.5, "mm"),
  column_gap = unit(0.5, "mm"),
  row_names_gp = gp_small,
  column_names_gp = gp_small,
  heatmap_legend_param = list(
    title_gp = gp_small,
    labels_gp = gp_small,
    legend_height = unit(20, "mm"),
    grid_height = unit(2, "mm"),
    grid_width = unit(4, "mm"),       
    direction = "horizontal"     
  )
)
w_mm <- 88
h_mm <- 110 

    # Save as PDF
pdf("/vol/projects/yzhang/500FG_aging/output/18_multi_omics_heatmap/multi_omics_heatmap_nature_gender.pdf", 
    width = w_mm / 25.4, height = h_mm / 25.4)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()
# Save as PNG
png("/vol/projects/yzhang/500FG_aging/output/18_multi_omics_heatmap/multi_omics_heatmap_nature_gender.png", 
        width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 300)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")

save.image("/vol/projects/yzhang/500FG_aging/output/18_multi_omics_heatmap/multi_omics_heatmap_nature_gender.Rdata")
