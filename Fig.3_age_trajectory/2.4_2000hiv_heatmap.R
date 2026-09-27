## Sup fig: 2000HIV validation cohort, same protein features as the 500FG
## panel (Fig 1A), fit to a single 80 mm column.
##


library(ComplexHeatmap)
library(dplyr)
library(circlize)
source("header.R")
COHORT_TARGET_WIDTH_MM <- 80   # this cohort's column width, not the shared TARGET_WIDTH_MM
TARGET_WIDTH_MM <- COHORT_TARGET_WIDTH_MM

## ---- paths (override via env vars; defaults are placeholders) -------------
cfg <- list(
  protein_file = Sys.getenv("HIV_PROTEIN_FILE", "data/raw/2000hiv_protein_trans.rds"),
  pheno_file   = Sys.getenv("HIV_PHENO_FILE", "data/raw/2000hiv_phenotype_trans.rds"),
  output_stem  = Sys.getenv("HIV_OUTPUT_STEM", "output/fig3c_heatmap_2000hiv")
)
dir.create(dirname(cfg$output_stem), recursive = TRUE, showWarnings = FALSE)

## Same protein feature list as the Fig 1A (500FG) panel; hardcoded here
## as in the source notebook rather than read from that panel's row-order
## RDS.
protein_500fg_names <- c(
  "IL18", "CCL25", "MCP_1", "IL_15RA", "CCL3", "IL6", "IL8", "CDCP1",
  "CST5", "Flt3L", "VEGFA", "MCP_2", "OPG", "CCL4", "HGF", "MCP_4",
  "CCL11", "MMP_1", "CCL28", "CXCL10", "CXCL11", "CXCL9", "CXCL1",
  "CXCL5", "CASP_8", "LAP_TGF_beta_1", "CD40", "IL_10RB", "CX3CL1",
  "TRANCE", "IL_12B", "TNFRSF9", "CD5", "TNFB", "CD8A", "SCF", "NT_3"
)

## ---- load, match 2000HIV column names to the 500FG names by prefix --------
protein <- readRDS(cfg$protein_file)
basicPhenos <- readRDS(cfg$pheno_file)

prefix_in_protein <- sub("_.*", "", colnames(protein))
prefix_500fg <- sub("_.*", "", protein_500fg_names)
full_names_matched <- colnames(protein)[match(prefix_500fg, prefix_in_protein)]
full_names_matched <- full_names_matched[!is.na(full_names_matched)]
protein_subset <- protein[, full_names_matched, drop = FALSE]

## ---- merge on common samples -----------------------------------------------
common_samples <- intersect(rownames(protein_subset), rownames(basicPhenos))
protein_filtered <- protein_subset[common_samples, , drop = FALSE]
basicPhenos_filtered <- as.data.frame(basicPhenos)[common_samples, , drop = FALSE]

merged_data <- data.frame(Sample_ID = common_samples, Age = basicPhenos_filtered$age, protein_filtered)

## ---- group by age (median per age), keep only complete rows ---------------
heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

valid_rows <- complete.cases(heatmap_data_sorted)
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted[valid_rows, , drop = FALSE])
age_sorted <- age_sorted[valid_rows]
heatmap_data_sorted$Age <- age_sorted

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
sample_counts <- age_counts$Count[match(age_sorted, age_counts$Age)]
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
age_col_fun <- colorRamp2(c(min(age_sorted), max(age_sorted)), c("#f2e8cf", "#a7c957"))

## ---- build + fit to the 80 mm column, then export --------------------------
build_heatmap_2000hiv <- function(p) {
  function() {
    protein_heatmap <- Heatmap(
      protein_data_zscore, name = "Proteomics Z-score",
      height = unit(nrow(protein_data_zscore) * p$ROW_HEIGHT_MM, "mm"),
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = TRUE, show_column_names = TRUE,
      row_names_gp = gpar(fontsize = p$FONT_PT), column_names_gp = gpar(fontsize = p$FONT_PT),
      col = protein_col_fun, row_dend_reorder = FALSE, width = unit(p$HEATMAP_WIDTH_MM, "mm"),
      top_annotation = HeatmapAnnotation(
        Age = anno_simple(age_sorted, col = age_col_fun, height = unit(p$AGE_BAR_HEIGHT_MM, "mm")),
        Age_Label = anno_text(age_labels, gp = gpar(fontsize = p$FONT_PT), rot = 0),
        annotation_name_side = "left", annotation_name_gp = gpar(fontsize = p$FONT_PT)
      ),
      left_annotation = rowAnnotation(
        Feature_Type = anno_block(gp = gpar(fill = "#457b9d", col = "white"), width = unit(p$LEFT_ANNO_WIDTH_MM, "mm"))
      ),
      rect_gp = gpar(col = NA), border = FALSE,
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"), labels_gp = gpar(fontsize = p$FONT_PT),
        direction = "horizontal", legend_width = unit(p$LEGEND_WIDTH_MM, "mm"), grid_height = unit(p$LEGEND_GRID_MM, "mm")
      )
    )

    feature_type_legend <- Legend(
      title = "Feature Type", title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
      at = c("Proteomics"), legend_gp = gpar(fill = c("#457b9d")),
      labels_gp = gpar(fontsize = p$FONT_PT), grid_height = unit(p$LEGEND_GRID_MM, "mm"), grid_width = unit(p$LEGEND_GRID_MM, "mm")
    )

    heatmap_combined <- make_sample_count_anno(sample_counts, p) %v% protein_heatmap

    draw(
      heatmap_combined, column_title = "2000HIV", column_title_side = "top",
      column_title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
      show_heatmap_legend = TRUE, show_annotation_legend = TRUE,
      annotation_legend_list = list(feature_type_legend),
      heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE
    )
  }
}

fit <- fit_to_journal_grid(base_params(), build_heatmap_2000hiv)
export_fitted_figure(cfg$output_stem, fit, caption_fn = NULL)
