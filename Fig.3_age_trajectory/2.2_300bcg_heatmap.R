## Sup Fig: 300BCG validation cohort, same features + row order as 500FG
##
## Requires header.R in the same directory, and the row-order
## RDS saved by heatmap_cytokine_protein.R (the 500FG panel this validates
## against) -- run that script first.
##

library(readxl)
source("header.R")

## ---- paths (override via env vars; defaults are placeholders) -------------
sup_cfg <- list(
  cytokine_file  = Sys.getenv("SUP_CYTOKINE_FILE", "data/raw/300bcg_cytokine.csv"),
  protein_file   = Sys.getenv("SUP_PROTEIN_FILE", "data/raw/300bcg_olink_baseline.RDS"),
  pheno_file     = Sys.getenv("SUP_PHENO_FILE", "data/300bcg_age_sex.xlsx"),
  row_order_file = Sys.getenv("HEATMAP_500FG_ROW_ORDER_FILE", "output/heatmap_500fg_cyto_pro_row_order.rds"),
  output_stem    = Sys.getenv("SUP_OUTPUT_STEM",
                               "output/sup_fig3b_300bcg_heatmap_zscore_samefeature")
)
dir.create(dirname(sup_cfg$output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- load, reorder protein columns to match the 500FG row order -----------
row_order_ref       <- readRDS(sup_cfg$row_order_file)
protein_500fg_order  <- row_order_ref$protein

cytokine_filtered <- read.csv(sup_cfg$cytokine_file, row.names = 1)  # all features kept, not FDR-filtered
protein_raw       <- readRDS(sup_cfg$protein_file)
pheno             <- read_excel(sup_cfg$pheno_file)

matching_columns <- protein_500fg_order[protein_500fg_order %in% colnames(protein_raw)]
protein_subset   <- protein_raw[, matching_columns, drop = FALSE]

## ---- merge on common samples, no imputation (rows with any NA are dropped
##      below instead -- matches the original analysis) ----------------------
## Note: the source script used an undefined `protein_filtered` here (the
## block that would have defined it was commented out, so it only worked by
## relying on a leftover variable from an interactively load()-ed .Rdata
## cache). `protein_subset`, computed just above, is what the visible code
## actually produces -- used here instead.
common_samples <- intersect(intersect(rownames(cytokine_filtered), rownames(protein_subset)),
                             pheno$PatientID)

cytokine_common <- cytokine_filtered[common_samples, , drop = FALSE]
protein_common  <- protein_subset[common_samples, , drop = FALSE]
pheno_common    <- as.data.frame(pheno[pheno$PatientID %in% common_samples, , drop = FALSE])
rownames(pheno_common) <- pheno_common$PatientID
pheno_common    <- pheno_common[common_samples, , drop = FALSE]

merged_data <- data.frame(Sample_ID = common_samples, Age = pheno_common$Age,
                           cytokine_common, protein_common)

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
  summarise(across(everything(), \(x) median(x, na.rm = TRUE))) %>%
  ungroup() %>%
  arrange(Age)

age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped_matrix <- as.matrix(dplyr::select(heatmap_data_grouped, -Age))

## Cytokine keeps its own column order; protein strictly follows the 500FG
## reference order so rows line up between the two cohort figures.
cytokine_cols_use <- colnames(cytokine_common)[colnames(cytokine_common) %in% colnames(heatmap_data_grouped_matrix)]
protein_cols_use  <- protein_500fg_order[protein_500fg_order %in% colnames(heatmap_data_grouped_matrix)]

cytokine_data <- t(heatmap_data_grouped_matrix[, cytokine_cols_use, drop = FALSE])
protein_data  <- t(heatmap_data_grouped_matrix[, protein_cols_use, drop = FALSE])

age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))
sample_counts  <- age_counts$Count[match(age_sorted, age_counts$Age)]
sample_counts[is.na(sample_counts)] <- 0

age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0,
  as.character(age_sorted), ""
)
age_labels[age_sorted == 70] <- ""  # hides an overlapping label in this cohort

## ---- z-score (0 for zero-variance rows, instead of NaN) -------------------
zscore_rows_safe <- function(m) {
  t(apply(m, 1, function(x) {
    x_sd <- sd(x, na.rm = TRUE)
    if (is.na(x_sd) || x_sd == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / x_sd
  }))
}
cytokine_data_zscore <- zscore_rows_safe(cytokine_data)
protein_data_zscore  <- zscore_rows_safe(protein_data)

cytokine_limit <- min(max(abs(range(cytokine_data_zscore, na.rm = TRUE))), 2.5)
protein_limit  <- min(max(abs(range(protein_data_zscore,  na.rm = TRUE))), 2.5)
cytokine_col_fun <- colorRamp2(seq(-cytokine_limit, cytokine_limit, length.out = 100),
                                colorRampPalette(c("#669bbc", "white", "#c1121f"))(100))
protein_col_fun  <- colorRamp2(seq(-protein_limit, protein_limit, length.out = 100),
                                colorRampPalette(c("#669bbc", "white", "#c1121f"))(100))
age_col_fun <- colorRamp2(c(min(age_sorted), max(age_sorted)), c("#f2e8cf", "#a7c957"))

stopifnot(
  length(age_sorted) == ncol(cytokine_data_zscore),
  length(age_sorted) == ncol(protein_data_zscore),
  length(sample_counts) == length(age_sorted)
)

## ---- build + export (fixed size -- matches the 500FG panel's proportions,
##      not run through fit_to_journal_grid()) --------------------------------
font_pt <- 5
row_height_mm <- 1.5
heatmap_width_mm <- 45

cytokine_heatmap <- Heatmap(
  cytokine_data_zscore, name = "Cytokine Expression",
  height = unit(nrow(cytokine_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = FALSE,
  row_names_gp = gpar(fontsize = font_pt), col = cytokine_col_fun,
  width = unit(heatmap_width_mm, "mm"),
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun, height = unit(2, "mm")),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = font_pt), rot = 0),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = font_pt)
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(gp = gpar(fill = "#e63946", col = "white"), width = unit(2, "mm"))
  ),
  rect_gp = gpar(col = NA), border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"), labels_gp = gpar(fontsize = font_pt),
    direction = "horizontal", legend_width = unit(15, "mm"), grid_height = unit(2, "mm")
  )
)

protein_heatmap <- Heatmap(
  protein_data_zscore, name = "Proteomics Z-score",
  height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt), column_names_gp = gpar(fontsize = font_pt),
  col = protein_col_fun, width = unit(heatmap_width_mm, "mm"),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(gp = gpar(fill = "#457b9d", col = "white"), width = unit(2, "mm"))
  ),
  rect_gp = gpar(col = NA), border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"), labels_gp = gpar(fontsize = font_pt),
    direction = "horizontal", legend_width = unit(15, "mm"), grid_height = unit(2, "mm")
  )
)

feature_type_legend_sup <- Legend(
  title = "Feature Type", title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Cytokine responses", "Proteomics"), legend_gp = gpar(fill = c("#e63946", "#457b9d")),
  labels_gp = gpar(fontsize = font_pt), grid_height = unit(2, "mm"), grid_width = unit(2, "mm")
)

sample_count_anno_sup <- HeatmapAnnotation(
  `Sample size` = anno_barplot(
    sample_counts, gp = gpar(fill = "#a8dadc", col = NA, lwd = 0), height = unit(0.8, "cm"),
    axis_param = list(at = c(0, max(sample_counts)), labels = c("0", max(sample_counts)),
                       gp = gpar(fontsize = font_pt)),
    border = FALSE
  ),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = font_pt)
)

heatmap_combined_sup <- sample_count_anno_sup %v% (cytokine_heatmap %v% protein_heatmap)

pdf(paste0(sup_cfg$output_stem, ".pdf"), width = 3.46, height = 3.9, onefile = FALSE)
draw(
  heatmap_combined_sup, column_title = "300BCG", column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE, show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_sup),
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE
)
dev.off()
