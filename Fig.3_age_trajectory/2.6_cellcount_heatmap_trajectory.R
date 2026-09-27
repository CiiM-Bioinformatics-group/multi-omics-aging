## Fig 1D+E (500FG): cell-count composition heatmap (3 immunotypes x 3
## feature clusters) and per-cluster age trajectory.
## Produces: heatmap_D_500fg_nature_final_legend.pdf / .png
##           trajectory_E_500fg_nature_final_legend.pdf / .png
##
## Requires header.R in the same directory -- it already
## defines CLUSTER_COLORS and common_ggplot_theme() for these two panels.
##

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
source("header.R")

## ---- paths (override via env vars; defaults are placeholders) -------------
cfg <- list(
  cellcount_sig_file     = Sys.getenv("CELLCOUNT_SIG_FILE", "data/cellcounts_age_rawless_gender_FDR.csv"),
  cellcount_raw_file     = Sys.getenv("CELLCOUNT_RAW_FILE", "data/raw/cellcounts_name_replace.csv"),
  cellcount_info_file    = Sys.getenv("CELLCOUNT_INFO_FILE", "data/500FG_cellcounts_info.txt"),
  pheno_file             = Sys.getenv("CELLCOUNT_PHENO_FILE", "data/Age_group_basicPhenos.csv"),
  heatmap_output_stem    = Sys.getenv("CELLCOUNT_HEATMAP_OUTPUT_STEM", "output/heatmap_D_500fg_nature_final_legend"),
  trajectory_output_stem = Sys.getenv("CELLCOUNT_TRAJECTORY_OUTPUT_STEM", "output/trajectory_E_500fg_nature_final_legend")
)
dir.create(dirname(cfg$heatmap_output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- load, filter to significant features, rename to display names --------
cellcount_sig <- read.csv(cfg$cellcount_sig_file, row.names = 1)
cellcount_sig <- cellcount_sig[cellcount_sig$sig == "sig", ]

cellcounts <- read.csv(cfg$cellcount_raw_file, row.names = 1, check.names = FALSE)
cellcounts_info <- read.table(cfg$cellcount_info_file)
colnames(cellcounts) <- cellcounts_info$finalName[match(colnames(cellcounts), cellcounts_info$trait)]

cellcount_exp <- cellcounts[, colnames(cellcounts) %in% cellcount_sig$feature]

pheno <- read.csv(cfg$pheno_file, row.names = 1)
rownames(pheno) <- pheno$ID_500fg

age_by_id <- setNames(pheno[rownames(cellcount_exp), "Age"], rownames(cellcount_exp))

## ---- cluster: columns (individuals) -> immunotype groups, rows (cell
## types) -> feature clusters -------------------------------------------------
mat <- t(as.matrix(cellcount_exp))
mat_scaled <- t(scale(t(mat)))

hc_col <- hclust(dist(t(mat_scaled), method = "euclidean"), method = "ward.D2")
immunotype <- cutree(hc_col, k = 3)
group_label <- factor(paste0("Group ", immunotype[colnames(mat)]), levels = paste0("Group ", 1:3))

hc_row <- hclust(dist(mat_scaled, method = "euclidean"), method = "ward.D2")
feature_cluster <- cutree(hc_row, k = 3)

## =============================================================================
## Panel D: composition heatmap, columns split by immunotype, rows split by
## feature cluster, Age annotated on top.
## =============================================================================
col_fun <- colorRamp2(
  seq(-5, 5, length.out = 9),
  rev(c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"))
)

top_anno <- HeatmapAnnotation(
  Age = anno_simple(age_by_id[colnames(mat)], col = colorRamp2(range(age_by_id, na.rm = TRUE), c("#f7f7f7", "#2166ac"))),
  annotation_name_gp = gpar(fontsize = 5)
)

ht <- Heatmap(
  mat_scaled, name = "Scaled cell proportion", top_annotation = top_anno,
  cluster_rows = hc_row, cluster_columns = hc_col,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = 3, column_split = 3,
  column_title = c("Group 1", "Group 2", "Group 3"), column_title_gp = gpar(fontsize = 5, fontface = "bold"),
  row_title = c("Cluster 1", "Cluster 2", "Cluster 3"), row_title_gp = gpar(fontsize = 5, fontface = "bold"),
  show_column_names = FALSE, show_row_names = TRUE,
  col = col_fun, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5),
    direction = "horizontal", legend_width = unit(1.5, "cm"), grid_height = unit(1, "mm")
  )
)

pdf(paste0(cfg$heatmap_output_stem, ".pdf"), width = 82.7 / 25.4, height = 103.4 / 25.4, onefile = FALSE)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

png(paste0(cfg$heatmap_output_stem, ".png"), width = 82.7 / 25.4, height = 103.4 / 25.4, units = "in", res = PNG_DPI)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

## =============================================================================
## Panel E: per-feature-cluster mean cell proportion vs. age (LOESS).
## =============================================================================
fc_vec_feature <- setNames(paste0("C", feature_cluster), rownames(mat))
mat_raw <- as.matrix(cellcount_exp)

## Proportions may already be in % (values up to ~100) or in [0, 1] --
## detect which, and scale to % only if needed.
already_pct <- max(mat_raw, na.rm = TRUE) > 1.5

traj_df <- mat_raw %>%
  as.data.frame() %>%
  tibble::rownames_to_column("individual") %>%
  pivot_longer(-individual, names_to = "cell_type", values_to = "proportion") %>%
  mutate(feature_cluster = fc_vec_feature[.data$cell_type]) %>%
  filter(!is.na(.data$feature_cluster)) %>%
  group_by(individual, feature_cluster) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE) * if (already_pct) 1 else 100, .groups = "drop") %>%
  mutate(Age = age_by_id[.data$individual]) %>%
  filter(is.finite(.data$Age)) %>%
  mutate(feature_cluster = factor(.data$feature_cluster, levels = c("C1", "C2", "C3")))

p_traj <- ggplot(traj_df, aes(x = Age, y = mean_prop, color = feature_cluster, fill = feature_cluster)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = CLUSTER_COLORS, name = "Feature cluster") +
  scale_fill_manual(values = CLUSTER_COLORS, name = "Feature cluster", guide = "none") +
  labs(x = "Age (years)", y = "Mean cell proportion (%)") +
  common_ggplot_theme(font_pt = 5) +
  theme(legend.position = "right")

ggsave(paste0(cfg$trajectory_output_stem, ".pdf"), p_traj, width = 64.9 / 25.4, height = 33.0 / 25.4)
ggsave(paste0(cfg$trajectory_output_stem, ".png"), p_traj, width = 64.9 / 25.4, height = 33.0 / 25.4,
       units = "in", dpi = PNG_DPI)
