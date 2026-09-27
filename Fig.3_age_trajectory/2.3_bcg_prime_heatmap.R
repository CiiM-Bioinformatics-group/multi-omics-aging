## Sup Fig: BCG-Prime cohort, protein features in the 500FG row order
## Produces: <PRIME_OUTPUT_STEM>.pdf
##
## Requires the row-order RDS saved by heatmap_cytokine_protein.R (the
## 500FG panel this validates against) -- run that script first.
##

library(dplyr)
library(ComplexHeatmap)
library(circlize)

## ---- paths (override via env vars; defaults are placeholders) -------------
cfg <- list(
  phenotype_file = Sys.getenv("PRIME_PHENOTYPE_FILE", "data/prime_phenotype_common.RDS"),
  olink_file     = Sys.getenv("PRIME_OLINK_FILE", "data/raw/prime_olink_raw_cleaned.RDS"),
  row_order_file = Sys.getenv("HEATMAP_500FG_ROW_ORDER_FILE", "output/heatmap_500fg_cyto_pro_row_order.rds"),
  output_stem    = Sys.getenv("PRIME_OUTPUT_STEM",
                               "output/sup_fig3c_heatmap_prime_zscore_same_feature"),
  out_rds        = Sys.getenv("PRIME_OUT_RDS", "output/heatmap_prime.rds")
)
dir.create(dirname(cfg$output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- load, reorder protein columns to match the 500FG row order -----------
phenotype <- readRDS(cfg$phenotype_file)

olink <- readRDS(cfg$olink_file)
olink_common <- olink[rownames(olink) %in% rownames(phenotype), ]
olink_common <- olink_common[match(rownames(phenotype), rownames(olink_common)), ]
olink_common_filled <- apply(olink_common, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

protein_500fg_order <- readRDS(cfg$row_order_file)$protein
matching_columns <- protein_500fg_order[protein_500fg_order %in% colnames(olink_common_filled)]
protein_subset <- olink_common_filled[, matching_columns, drop = FALSE]

merged_data <- data.frame(Sample_ID = rownames(phenotype), Age = phenotype$Age, protein_subset)

## ---- group by age (median per age), keep only complete rows ---------------
heatmap_data   <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]
sorted_idx     <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_idx, , drop = FALSE]
age_sorted_raw <- merged_data$Age[sorted_idx]

valid_rows <- complete.cases(heatmap_data_sorted)
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted[valid_rows, , drop = FALSE])
age_sorted_raw <- age_sorted_raw[valid_rows]
heatmap_data_sorted$Age <- age_sorted_raw

heatmap_data_grouped <- heatmap_data_sorted %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)

age_sorted <- heatmap_data_grouped$Age
protein_data <- t(as.matrix(select(heatmap_data_grouped, -Age)))

age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))
sample_counts  <- age_counts$Count[match(age_sorted, age_counts$Age)]
sample_counts[is.na(sample_counts)] <- 0
stopifnot(length(sample_counts) == length(age_sorted))

age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0,
  as.character(age_sorted), ""
)

## ---- z-score (0 for zero-variance rows, instead of NaN) -------------------
zscore_data <- function(data_matrix) {
  t(apply(data_matrix, 1, function(x) {
    x_sd <- sd(x, na.rm = TRUE)
    if (is.na(x_sd) || x_sd == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / x_sd
  }))
}
protein_data_zscore <- zscore_data(protein_data)

protein_limit <- min(max(abs(range(protein_data_zscore, na.rm = TRUE))), 2.5)
protein_col_fun <- colorRamp2(seq(-protein_limit, protein_limit, length.out = 100),
                               colorRampPalette(c("#669bbc", "white", "#c1121f"))(100))
age_col_fun <- colorRamp2(c(min(age_sorted), max(age_sorted)), c("#cfd796", "#a7c957"))

## ---- build + export (fixed size -- matches the source figure's proportions) ----
font_pt <- 5
row_height_mm <- 1.5
heatmap_width_mm <- 45

protein_heatmap <- Heatmap(
  protein_data_zscore, name = "Proteomics Z-score",
  height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt), column_names_gp = gpar(fontsize = font_pt),
  col = protein_col_fun, width = unit(heatmap_width_mm, "mm"),
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun, height = unit(2, "mm")),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = font_pt), rot = 0),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = font_pt)
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(gp = gpar(fill = "#457b9d", col = "white"), width = unit(2, "mm"))
  ),
  rect_gp = gpar(col = NA), border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"), labels_gp = gpar(fontsize = font_pt),
    direction = "horizontal", legend_width = unit(15, "mm"), grid_height = unit(2, "mm")
  )
)

feature_type_legend <- Legend(
  title = "Feature Type", title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Proteomics"), legend_gp = gpar(fill = c("#457b9d")),
  labels_gp = gpar(fontsize = font_pt), grid_height = unit(2, "mm"), grid_width = unit(2, "mm")
)

sample_count_anno <- HeatmapAnnotation(
  `Sample size` = anno_barplot(
    sample_counts, gp = gpar(fill = "#a8dadc", col = NA, lwd = 0), height = unit(0.8, "cm"),
    axis_param = list(at = c(0, max(sample_counts)), labels = c("0", max(sample_counts)),
                       gp = gpar(fontsize = font_pt)),
    border = FALSE
  ),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = font_pt)
)

heatmap_combined <- sample_count_anno %v% protein_heatmap

pdf(paste0(cfg$output_stem, ".pdf"), width = 3.46, height = 3.5, onefile = FALSE)
draw(
  heatmap_combined, column_title = "BCG-Prime", column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE, show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend),
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE
)
dev.off()

saveRDS(
  list(protein_data = protein_data, sample_counts = sample_counts,
       age_sorted = age_sorted, age_labels = age_labels, age_col_fun = age_col_fun),
  cfg$out_rds
)
