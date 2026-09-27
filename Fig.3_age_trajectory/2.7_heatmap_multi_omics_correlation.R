## Multi-omics correlation heatmap (500FG): feature-feature Spearman
## correlation across age-associated features from every layer.
## Produces: multi_omics_heatmap_nature_gender.pdf / .png

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

## ---- paths (override via env vars; defaults are placeholders) -------------
sig <- list(
  methylation     = Sys.getenv("CORR_METHYLATION_SIG_FILE", "data/methylation_age_liner_gender_FDR.RDS"),
  olink           = Sys.getenv("CORR_OLINK_SIG_FILE", "data/olink_age_liner_gender_FDR.csv"),
  cytokine        = Sys.getenv("CORR_CYTOKINE_SIG_FILE", "data/cytokine_age_liner_gender_FDR.csv"),
  metabolite      = Sys.getenv("CORR_METABOLITE_SIG_FILE", "data/metabolite_age_linear_gender_FDR.csv"),
  hormone         = Sys.getenv("CORR_HORMONE_SIG_FILE", "data/hormone_age_liner_gender_FDR.csv"),
  cellcounts      = Sys.getenv("CORR_CELLCOUNTS_SIG_FILE", "data/cellcounts_age_rawless_gender_FDR.csv"),
  immunoglobulin  = Sys.getenv("CORR_IMMUNOGLOBULIN_SIG_FILE", "data/immunoglobulin_age_gender_FDR.csv"),
  microbiome      = Sys.getenv("CORR_MICROBIOME_SIG_FILE", "data/microbiome_age_gender_FDR.csv")
)
raw <- list(
  methylation     = Sys.getenv("CORR_METHYLATION_RAW_FILE", "data/raw/Mvalue_trans.rds"),
  olink           = Sys.getenv("CORR_OLINK_RAW_FILE", "data/raw/olink_NPX.RDS"),
  cytokine        = Sys.getenv("CORR_CYTOKINE_RAW_FILE", "data/raw/cytokine_filled_trans.csv"),
  metabolite      = Sys.getenv("CORR_METABOLITE_RAW_FILE", "data/raw/metabolite_trans.rds"),
  hormone         = Sys.getenv("CORR_HORMONE_RAW_FILE", "data/raw/hormone_log2_levels.txt"),
  cellcounts      = Sys.getenv("CORR_CELLCOUNTS_RAW_FILE", "data/raw/cellcounts_name_replace.csv"),
  microbiome      = Sys.getenv("CORR_MICROBIOME_RAW_FILE", "data/raw/microbiome_pathways.txt"),
  immunoglobulin  = Sys.getenv("CORR_IMMUNOGLOBULIN_RAW_FILE", "data/raw/immunoglobulin_log2_levels.txt")
)
output_stem <- Sys.getenv("CORR_OUTPUT_STEM", "output/sup_fig4_multi_omics_heatmap_nature_gender")
dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- 1. significant features per layer -------------------------------------
filter_sig <- function(df, layer, is_filtered = FALSE) {
  if (!is_filtered) df <- df[df$sig == "sig", ]
  data.frame(feature = df$feature, estimate = df$estimate, layer = layer)
}

methylation <- readRDS(sig$methylation) %>% rename(feature = Mvalue)
methylation_filtered <- methylation[methylation$padj < 0.0001 & abs(methylation$estimate) > 0.02, ]

olink          <- read.csv(sig$olink, row.names = 1)
cytokine       <- read.csv(sig$cytokine, row.names = 1)
metabolite     <- read.csv(sig$metabolite, row.names = 1)
hormone        <- read.csv(sig$hormone, row.names = 1)
cellcounts     <- read.csv(sig$cellcounts, row.names = 1)
immunoglobulin <- read.csv(sig$immunoglobulin, row.names = 1)
microbiome     <- read.csv(sig$microbiome, row.names = 1)

feature_info <- rbind(
  filter_sig(methylation_filtered, "DNA methylation", TRUE),
  filter_sig(olink, "Proteomics"),
  filter_sig(cytokine, "Cytokine response"),
  filter_sig(metabolite, "Metabolomics"),
  filter_sig(cellcounts, "Immune cell counts"),
  filter_sig(microbiome, "Microbiomes"),
  filter_sig(hormone, "Circulating endocrine traits"),
  filter_sig(immunoglobulin, "Circulating immune traits")
)

## ---- 2. subset each layer's raw data to its significant features ----------
select_features <- function(raw_df, sig_features) {
  raw_df[, colnames(raw_df) %in% sig_features, drop = FALSE]
}

## Methylation's raw matrix (854k CpGs) is only ever needed for the already
## -filtered feature set, so it's dropped from memory right after subsetting.
Mvalue_raw <- readRDS(raw$methylation)
Mvalue_sub <- Mvalue_raw[, colnames(Mvalue_raw) %in% methylation_filtered$feature]
rm(Mvalue_raw)

olink_raw          <- readRDS(raw$olink)
cytokine_raw       <- read.csv(raw$cytokine, row.names = 1)
metabolite_raw     <- readRDS(raw$metabolite)
hormone_raw        <- read.table(raw$hormone, header = TRUE, sep = "", stringsAsFactors = FALSE)
cellcounts_raw     <- read.csv(raw$cellcounts, row.names = 1, check.names = FALSE)
microbiome_raw     <- read.table(raw$microbiome, header = TRUE, sep = "", stringsAsFactors = FALSE)
immunoglobulin_raw <- read.table(raw$immunoglobulin, header = TRUE, sep = "", stringsAsFactors = FALSE)

olink_sub          <- select_features(olink_raw, feature_info$feature[feature_info$layer == "Proteomics"])
cytokine_sub       <- select_features(cytokine_raw, feature_info$feature[feature_info$layer == "Cytokine response"])
metabolite_sub     <- select_features(metabolite_raw, feature_info$feature[feature_info$layer == "Metabolomics"])
cellcounts_sub     <- select_features(cellcounts_raw, feature_info$feature[feature_info$layer == "Immune cell counts"])
microbiome_sub     <- select_features(microbiome_raw, feature_info$feature[feature_info$layer == "Microbiomes"])
hormone_sub        <- select_features(hormone_raw, feature_info$feature[feature_info$layer == "Circulating endocrine traits"])
immunoglobulin_sub <- select_features(immunoglobulin_raw, feature_info$feature[feature_info$layer == "Circulating immune traits"])

## ---- 3-4. common samples, merge into one feature x sample matrix ----------
common_samples <- Reduce(intersect, list(
  rownames(Mvalue_sub), rownames(olink_sub), rownames(cytokine_sub), rownames(metabolite_sub),
  rownames(cellcounts_sub), rownames(microbiome_sub), rownames(hormone_sub), rownames(immunoglobulin_sub)
))

raw_data <- cbind(
  Mvalue_sub[common_samples, ], olink_sub[common_samples, ], cytokine_sub[common_samples, ],
  metabolite_sub[common_samples, ], cellcounts_sub[common_samples, ], microbiome_sub[common_samples, ],
  hormone_sub[common_samples, ], immunoglobulin_sub[common_samples, ]
)

feature_info <- feature_info[feature_info$feature %in% colnames(raw_data), ]
layer_order <- c("DNA methylation", "Proteomics", "Cytokine response", "Metabolomics",
                  "Immune cell counts", "Microbiomes", "Circulating endocrine traits", "Circulating immune traits")
feature_info$layer <- factor(feature_info$layer, levels = layer_order)

## ---- 5. feature-feature Spearman correlation, FDR-adjusted p-values -------
cor_matrix <- cor(raw_data, use = "pairwise.complete.obs", method = "spearman")

layer_annotation <- feature_info$layer
names(layer_annotation) <- feature_info$feature
age_cor <- ifelse(feature_info$estimate > 0, "Positive", "Negative")
names(age_cor) <- feature_info$feature

feature_order <- feature_info$feature[order(feature_info$layer)]  # blocks by layer; clustered within each below

## p-value via t-approximation for Spearman, then BH-adjusted over unique pairs
n <- nrow(raw_data)
t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
t_stat[abs(cor_matrix) >= 1 | is.na(cor_matrix)] <- 0
p_matrix <- 2 * pt(-abs(t_stat), n - 2)
p_matrix[is.na(cor_matrix)] <- NA

p_adj_matrix <- matrix(NA, nrow(p_matrix), ncol(p_matrix))
p_adj_vec <- p.adjust(p_matrix[upper.tri(p_matrix)], method = "BH")
p_adj_matrix[upper.tri(p_adj_matrix)] <- p_adj_vec
p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
diag(p_adj_matrix) <- NA
p_adj_matrix[is.na(p_matrix)] <- NA

## Heatmap value: correlation direction, strength from FDR -> signed -log10(FDR),
## capped so a handful of near-zero FDRs don't blow out the color scale.
signed_logp <- sign(cor_matrix) * (-log10(p_adj_matrix))
signed_logp[is.infinite(signed_logp)] <- NA
cap <- 10  # -log10(1e-10)
signed_logp[signed_logp >  cap] <-  cap
signed_logp[signed_logp < -cap] <- -cap
signed_logp_ordered <- signed_logp[feature_order, feature_order]

## ---- 6. build + export ------------------------------------------------------
## Same layer palette as the other age-explained-variance figures, so a
## given layer reads as the same color everywhere in the manuscript.
layer_colors <- c(
  "DNA methylation"               = "#8DD3C7",
  "Proteomics"                    = "#FB8072",
  "Cytokine response"             = "#80B1D3",
  "Metabolomics"                  = "#E9C46A",
  "Immune cell counts"            = "#BEBADA",
  "Microbiomes"                   = "#FCCDE5",
  "Circulating endocrine traits"  = "#D98C6B",
  "Circulating immune traits"     = "#5FAF9D"
)
age_colors <- c("Positive" = "#FB8072", "Negative" = "#80B1D3")
col_fun <- colorRamp2(c(-1, 0, 1), c("#669bbc", "white", "#c1121f"))

gp_small    <- gpar(fontsize = 5)
grid_small  <- unit(2, "mm")
anno_legend_param <- list(
  `Age correlation` = list(title_gp = gp_small, labels_gp = gp_small,
                            grid_height = grid_small, grid_width = grid_small, nrow = 2),
  Layer             = list(title_gp = gp_small, labels_gp = gp_small,
                            grid_height = grid_small, grid_width = grid_small, nrow = 4)
)

row_ha <- rowAnnotation(
  `Age correlation` = age_cor[feature_order], Layer = layer_annotation[feature_order],
  col = list(`Age correlation` = age_colors, Layer = layer_colors),
  show_legend = FALSE, show_annotation_name = FALSE,
  annotation_name_gp = gp_small, simple_anno_size = unit(2, "mm")
)
col_ha <- HeatmapAnnotation(
  `Age correlation` = age_cor[feature_order], Layer = layer_annotation[feature_order],
  col = list(`Age correlation` = age_colors, Layer = layer_colors),
  show_legend = TRUE, show_annotation_name = FALSE,
  annotation_name_gp = gp_small, annotation_legend_param = anno_legend_param,
  simple_anno_size = unit(2, "mm")
)

ht <- Heatmap(
  signed_logp_ordered, name = "-log10 FDR", col = col_fun,
  top_annotation = col_ha, left_annotation = row_ha,
  show_row_names = FALSE, show_column_names = FALSE,
  cluster_rows = TRUE, cluster_columns = TRUE,
  cluster_row_slices = FALSE, cluster_column_slices = FALSE,
  row_split = layer_annotation[feature_order], column_split = layer_annotation[feature_order],
  row_title = NULL, column_title = NULL, border = FALSE,
  row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"),
  row_names_gp = gp_small, column_names_gp = gp_small,
  heatmap_legend_param = list(
    title_gp = gp_small, labels_gp = gp_small, legend_height = unit(20, "mm"),
    grid_height = unit(2, "mm"), grid_width = unit(4, "mm"), direction = "horizontal"
  )
)

w_mm <- 88
h_mm <- 110

pdf(paste0(output_stem, ".pdf"), width = w_mm / 25.4, height = h_mm / 25.4)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

png(paste0(output_stem, ".png"), width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
