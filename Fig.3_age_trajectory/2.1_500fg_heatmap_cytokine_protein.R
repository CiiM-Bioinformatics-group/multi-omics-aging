## Fig 1A heatmap: cytokine response + proteomics, age-median z-score (500FG)
## Produces: heatmap_500fg_nature_final_legend.pdf / .png
##
## Requires header.R in the same directory 
##

library(dplyr)
source("header.R")

## ---- paths (override via env vars; defaults are placeholders) -------------
cfg <- list(
  cytokine_file     = Sys.getenv("HEATMAP_CYTOKINE_FILE", "data/raw/cytokine_filled_trans.csv"),
  protein_file      = Sys.getenv("HEATMAP_PROTEIN_FILE", "data/raw/olink_NPX.RDS"),
  pheno_file        = Sys.getenv("HEATMAP_PHENO_FILE", "data/Age_group_basicPhenos.csv"),
  cytokine_sig_file = Sys.getenv("HEATMAP_CYTOKINE_SIG_FILE", "data/cytokine_age_liner_gender_FDR.csv"),
  protein_sig_file  = Sys.getenv("HEATMAP_PROTEIN_SIG_FILE", "data/olink_age_liner_gender_FDR.csv"),
  output_stem       = Sys.getenv("HEATMAP_OUTPUT_STEM", "output/fig3a_heatmap_500fg_final_legend")
)
dir.create(dirname(cfg$output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- load, keep only sex-adjusted FDR-significant features ----------------
cytokine    <- read.csv(cfg$cytokine_file, row.names = 1)
protein     <- as.data.frame(readRDS(cfg$protein_file))
basicPhenos <- read.csv(cfg$pheno_file, row.names = 1)

cytokine_sig <- read.csv(cfg$cytokine_sig_file, row.names = 1)
protein_sig  <- read.csv(cfg$protein_sig_file, row.names = 1)

cytokine_features <- cytokine_sig$feature[cytokine_sig$sig == "sig"]
protein_features  <- protein_sig$feature[protein_sig$sig == "sig"]
cytokine_filtered <- cytokine[, colnames(cytokine) %in% cytokine_features, drop = FALSE]
protein_filtered  <- protein[, colnames(protein) %in% protein_features, drop = FALSE]

## ---- merge on common samples, mean-impute remaining NAs -------------------
common_samples <- intersect(intersect(rownames(cytokine_filtered), rownames(protein_filtered)),
                             basicPhenos$ID_500fg)

cytokine_common    <- cytokine_filtered[common_samples, , drop = FALSE]
protein_common     <- protein_filtered[common_samples, , drop = FALSE]
basicPhenos_common <- basicPhenos[match(common_samples, basicPhenos$ID_500fg), ]

merged_data <- data.frame(Sample_ID = common_samples, Age = basicPhenos_common$Age,
                           cytokine_common, protein_common)
merged_data <- merged_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

## ---- group samples by age (median per age), z-score each feature ----------
heatmap_data   <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]
sorted_idx     <- order(merged_data$Age, na.last = FALSE)
age_sorted_raw <- merged_data$Age[sorted_idx]

heatmap_data_grouped <- heatmap_data[sorted_idx, , drop = FALSE] %>%
  mutate(Age = age_sorted_raw) %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)

age_sorted               <- heatmap_data_grouped$Age
heatmap_data_transposed  <- t(as.matrix(select(heatmap_data_grouped, -Age)))
cytokine_data <- heatmap_data_transposed[colnames(cytokine_common), , drop = FALSE]
protein_data  <- heatmap_data_transposed[colnames(protein_common), , drop = FALSE]

zscore_rows <- function(m) t(apply(m, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
cytokine_data_zscore <- zscore_rows(cytokine_data)
protein_data_zscore  <- zscore_rows(protein_data)

age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))
sample_counts  <- age_counts$Count[match(age_sorted, age_counts$Age)]
sample_counts[is.na(sample_counts)] <- 0

## ---- cytokine row metadata, shown as color strips instead of row names ----
n_500fg <- length(common_samples)
age_labels_500fg <- make_age_labels(age_sorted)

SHOW_CYTOKINE_ROW_NAMES <- FALSE  # cytokine rows are labelled via the strips below

cytokine_labels_clean <- clean_cytokine_label(rownames(cytokine_data_zscore))
protein_labels_clean  <- gsub("_", "-", rownames(protein_data_zscore))

## Raw feature names follow cytokine_stimulation_cellsystem_duration, e.g.
## "IL6_C.albicanshyphae_PBMC_24h" -- split on the raw name, not the
## cleaned display string.
meta_cyt <- tidyr::separate(
  data.frame(feature = rownames(cytokine_data_zscore), stringsAsFactors = FALSE),
  feature, into = c("cytokine", "stimulation", "cell_system", "duration"),
  sep = "_", remove = FALSE, fill = "right", extra = "merge"
)

cytokine_display_map <- c(IL6 = "IL-6", IL17 = "IL-17", IL22 = "IL-22")
meta_cyt$cytokine_disp <- ifelse(meta_cyt$cytokine %in% names(cytokine_display_map),
                                  cytokine_display_map[meta_cyt$cytokine], meta_cyt$cytokine)

## Stimulation has more raw prep-level categories than organisms; grouped
## by organism for a readable legend.
stimulation_group_map <- c(
  "C.conidiaHK" = "C. albicans", "C.albicansconidia" = "C. albicans", "C.albicanshyphae" = "C. albicans",
  "B.burgdorferi" = "Borrelia", "Borreliamix" = "Borrelia",
  "A.fumigatusconidiaSerum" = "A. fumigatus", "Bacteroides" = "Bacteroides",
  "MTB" = "MTB", "PHA" = "PHA", "S.aureus" = "S. aureus"
)
meta_cyt$stimulation_group <- unname(stimulation_group_map[meta_cyt$stimulation])

duration_to_hours <- function(x) {
  hrs <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  ifelse(grepl("day", x, ignore.case = TRUE), hrs * 24, hrs)
}

cytokine_levels    <- unique(meta_cyt$cytokine_disp)
cellsystem_levels  <- unique(meta_cyt$cell_system)
stimulation_levels <- unique(meta_cyt$stimulation_group)
duration_levels    <- unique(meta_cyt$duration)
duration_levels    <- duration_levels[order(duration_to_hours(duration_levels))]
meta_cyt$cytokine_disp <- factor(meta_cyt$cytokine_disp, levels = cytokine_levels)

## Muted, warm-to-cool palettes, each ramped light -> dark within its own track.
col_cytokine    <- setNames(colorRampPalette(c("#E5D8BE", "#AC8B52"))(length(cytokine_levels)), cytokine_levels)
col_cellsystem  <- setNames(colorRampPalette(c("#CFDDD2", "#668773"))(length(cellsystem_levels)), cellsystem_levels)
col_stimulation <- setNames(colorRampPalette(c("#D5E1E6", "#B9CCD5", "#9DB6C3", "#7899AC",
                                                "#5D8196", "#496778", "#354D5C"))(length(stimulation_levels)),
                             stimulation_levels)
col_duration    <- setNames(colorRampPalette(c("#DDD1DF", "#927397"))(length(duration_levels)), duration_levels)

## ---- build the panel; called twice by fit_to_journal_grid() (measure + fit) ----
build_heatmap_A <- function(p) {
  cytokine_labels_disp <- render_labels(cytokine_labels_clean, fontsize = p$FONT_PT)
  protein_labels_disp  <- render_labels(protein_labels_clean,  fontsize = p$FONT_PT)

  cytokine_left_anno <- rowAnnotation(
    `Feature type` = anno_block(gp = gpar(fill = "#e63946", col = "white"),
                                 width = unit(p$LEFT_ANNO_WIDTH_MM, "mm"))
  )
  cytokine_right_anno <- rowAnnotation(
    Cytokine = meta_cyt$cytokine_disp, `Cell system` = meta_cyt$cell_system,
    Stimulation = meta_cyt$stimulation_group, Duration = meta_cyt$duration,
    col = list(Cytokine = col_cytokine, `Cell system` = col_cellsystem,
               Stimulation = col_stimulation, Duration = col_duration),
    simple_anno_size = unit(p$META_ANNO_WIDTH_MM, "mm"),
    gap = unit(p$META_ANNO_GAP_MM, "mm"),
    show_legend = FALSE,  # built as separate Legend() objects below
    show_annotation_name = TRUE, annotation_name_side = "top",
    annotation_name_rot = 45, annotation_name_gp = gpar(fontsize = p$FONT_PT)
  )
  ## Reused by protein_heatmap's right_annotation so the two heatmap bodies
  ## stay aligned when stacked with %v%.
  total_right_width_mm <- 4 * p$META_ANNO_WIDTH_MM + 3 * p$META_ANNO_GAP_MM

  cytokine_heatmap <- Heatmap(
    cytokine_data_zscore,
    name = "cytokine_zscore",  # must differ from protein_heatmap's name
    col = zscore_col_fun,
    cluster_rows = TRUE, cluster_columns = FALSE,
    row_dend_width = unit(p$DEND_WIDTH_MM, "mm"),
    show_row_names = SHOW_CYTOKINE_ROW_NAMES, row_labels = cytokine_labels_disp,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = p$FONT_PT), column_names_gp = gpar(fontsize = p$FONT_PT),
    width = unit(p$HEATMAP_WIDTH_MM, "mm"),
    height = unit(nrow(cytokine_data_zscore) * p$ROW_HEIGHT_MM, "mm"),
    top_annotation = HeatmapAnnotation(
      Age = anno_simple(age_sorted, col = age_col_fun, height = unit(p$AGE_BAR_HEIGHT_MM, "mm")),
      Age_Label = anno_text(age_labels_500fg, gp = gpar(fontsize = p$FONT_PT), rot = 0),
      annotation_name_side = "left", annotation_name_gp = gpar(fontsize = p$FONT_PT)
    ),
    left_annotation = cytokine_left_anno, right_annotation = cytokine_right_anno,
    rect_gp = gpar(col = NA), border = FALSE,
    heatmap_legend_param = list(
      title = ZSCORE_LEGEND_NAME, title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
      labels_gp = gpar(fontsize = p$FONT_PT), direction = "horizontal",
      legend_width = unit(p$LEGEND_WIDTH_MM, "mm"), grid_height = unit(p$LEGEND_GRID_MM, "mm"),
      at = c(-Z_LIMIT, 0, Z_LIMIT)
    )
  )

  protein_heatmap <- Heatmap(
    protein_data_zscore,
    name = "protein_zscore",
    col = zscore_col_fun,
    cluster_rows = TRUE, cluster_columns = FALSE, row_dend_reorder = FALSE,
    row_dend_width = unit(p$DEND_WIDTH_MM, "mm"),
    show_row_names = FALSE, row_labels = protein_labels_disp,  # labelled via right_annotation instead
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = p$FONT_PT), column_names_gp = gpar(fontsize = p$FONT_PT),
    width = unit(p$HEATMAP_WIDTH_MM, "mm"),
    height = unit(nrow(protein_data_zscore) * p$ROW_HEIGHT_MM, "mm"),
    left_annotation = rowAnnotation(
      `Feature type` = anno_block(gp = gpar(fill = "#457b9d", col = "white"),
                                   width = unit(p$LEFT_ANNO_WIDTH_MM, "mm"))
    ),
    right_annotation = rowAnnotation(
      Feature = anno_text(
        protein_labels_disp, gp = gpar(fontsize = p$FONT_PT),
        width = max(unit(total_right_width_mm, "mm"),
                    ComplexHeatmap::max_text_width(protein_labels_disp, gp = gpar(fontsize = p$FONT_PT)) + unit(2, "mm"))
      ),
      show_annotation_name = FALSE
    ),
    rect_gp = gpar(col = NA), border = FALSE,
    show_heatmap_legend = FALSE  # cytokine_heatmap's legend is reused
  )

  feature_type_legend <- Legend(
    title = "Feature type", title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
    at = c("Cytokine response", "Proteomics"), legend_gp = gpar(fill = c("#e63946", "#457b9d")),
    labels_gp = gpar(fontsize = p$FONT_PT), grid_width = unit(p$LEGEND_GRID_MM, "mm"), ncol = 2
  )
  zscore_legend <- Legend(
    col_fun = zscore_col_fun, title = ZSCORE_LEGEND_NAME, at = c(-Z_LIMIT, 0, Z_LIMIT),
    direction = "horizontal", title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
    labels_gp = gpar(fontsize = p$FONT_PT), legend_width = unit(p$LEGEND_WIDTH_MM, "mm"),
    grid_height = unit(p$LEGEND_GRID_MM, "mm")
  )

  ## Each strip's own legend, entries flowing left-to-right (max 4/row).
  horiz_ncol <- function(n) min(n, 4)
  legend_common <- list(
    title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"), labels_gp = gpar(fontsize = p$FONT_PT),
    grid_width = unit(p$LEGEND_GRID_MM, "mm"), grid_height = unit(p$LEGEND_GRID_MM, "mm"), by_row = TRUE
  )
  cytokine_legend    <- do.call(Legend, c(list(title = "Cytokine", at = names(col_cytokine),
                                 legend_gp = gpar(fill = col_cytokine), ncol = horiz_ncol(length(col_cytokine))), legend_common))
  cellsystem_legend  <- do.call(Legend, c(list(title = "Cell system", at = names(col_cellsystem),
                                 legend_gp = gpar(fill = col_cellsystem), ncol = horiz_ncol(length(col_cellsystem))), legend_common))
  stimulation_legend <- do.call(Legend, c(list(title = "Stimulation", at = names(col_stimulation),
                                 legend_gp = gpar(fill = col_stimulation), ncol = horiz_ncol(length(col_stimulation))), legend_common))
  duration_legend    <- do.call(Legend, c(list(title = "Duration", at = names(col_duration),
                                 legend_gp = gpar(fill = col_duration), ncol = horiz_ncol(length(col_duration))), legend_common))

  ## Capped at TARGET_WIDTH_MM; wraps to a new row instead of overflowing.
  legend_pack <- packLegend(
    feature_type_legend, zscore_legend, cytokine_legend, cellsystem_legend, stimulation_legend, duration_legend,
    direction = "horizontal", max_width = unit(TARGET_WIDTH_MM - 2, "mm"),
    gap = unit(2, "mm"), row_gap = unit(1.5, "mm")
  )

  sample_count_anno <- make_sample_count_anno(sample_counts, p)
  heatmap_A <- sample_count_anno %v% (cytokine_heatmap %v% protein_heatmap)

  function() {
    draw(
      heatmap_A, column_title = "500FG", column_title_side = "top",
      column_title_gp = gpar(fontsize = p$FONT_PT, fontface = "bold"),
      show_heatmap_legend = FALSE,  # replaced by zscore_legend in legend_pack
      show_annotation_legend = TRUE, annotation_legend_list = list(legend_pack),
      annotation_legend_side = "bottom", ht_gap = unit(c(p$HT_GAP_MM, p$HT_GAP_MM / 4), "mm")
    )
  }
}

## ---- fit to the journal grid and export ------------------------------------
p_base <- base_params()
p_base$META_ANNO_WIDTH_MM <- 1.8
p_base$META_ANNO_GAP_MM   <- 0.3

fit_A <- fit_to_journal_grid(p_base, build_heatmap_A, height_scale = 2 / 3)
export_fitted_figure(cfg$output_stem, fit_A, caption_fn = NULL)
