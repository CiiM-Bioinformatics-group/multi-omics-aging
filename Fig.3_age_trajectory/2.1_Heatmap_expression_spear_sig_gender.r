library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)  
library(dendextend) 
library(NbClust)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(viridis)
library(readxl)

load("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_bar_zscore_spears_sig_feature_gender_nature.Rdata")

cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/input/cytokine/filtered_cytokine_filled_trans.csv",row.names=1)
head(cytokine)
dim(cytokine)

protein <- readRDS("/vol/projects/CIIM/cohorts_old/500FG/olinkData/NPX_FG500.RDS")
head(protein)
dim(protein)

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv",row.names=1)
head(basicPhenos)
dim(basicPhenos)

cytokine_sig <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv",row.names=1)
head(cytokine_sig)
dim(cytokine_sig)

protein_sig <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv",row.names=1)
head(protein_sig)
dim(protein_sig)

cytokine_sig <- cytokine_sig[cytokine_sig$sig=="sig",]
cytokine_features <- cytokine_sig$feature
cytokine_filtered <- cytokine[, colnames(cytokine) %in% cytokine_features]
head(cytokine_filtered)
dim(cytokine_filtered)

protein_sig <- protein_sig[protein_sig$sig == "sig", ]
protein_features <- protein_sig$feature
protein_filtered <- protein[, colnames(protein) %in% protein_features]
head(protein_filtered)
dim(protein_filtered)

common_samples <- intersect(intersect(rownames(cytokine_filtered), 
                                      rownames(protein_filtered)), 
                            basicPhenos$ID_500fg)

# Filter out common samples
cytokine_common <- cytokine_filtered[common_samples, , drop = FALSE]
protein_common <- protein_filtered[common_samples, , drop = FALSE]
basicPhenos_common <- basicPhenos[basicPhenos$ID_500fg %in% common_samples, , drop = FALSE]

# Reset row names of basicPhenos_filtered to match common_samples
rownames(basicPhenos_common) <- basicPhenos_common$ID_500fg

# Ensure basicPhenos_filtered only contains common samples
basicPhenos_common <- basicPhenos_common[common_samples, , drop = FALSE]

# Merge data (retain all features)
merged_data <- data.frame(
  Sample_ID = common_samples,
  Age = basicPhenos_common$Age,  # Age alignment
  cytokine_common,  # All cytokine features
  protein_common    # All protein features
)

# View merged data
head(merged_data)
dim(merged_data)

write.csv(merged_data,"/vol/projects/yzhang/500FG_aging/output/14_AA/500fg_heatmap_merge_data_gender.csv",row.names=TRUE)

protein_cyto_sig_feature <- colnames(merged_data)
protein_cyto_sig_feature <- protein_cyto_sig_feature[-c(1, 2)]
head(protein_cyto_sig_feature)
length(protein_cyto_sig_feature)
#write.csv(protein_cyto_sig_feature,"/vol/projects/yzhang/500FG_aging/output/14_AA/protein_cyto_sig_feature.csv")

##deal with missing value
sum(is.na(merged_data))
merged_data <- merged_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
sum(is.na(merged_data))

##nature style
# Prepare heatmap data (exclude Sample_ID and Age columns)
heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

# Sort data by age
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

# Group by age and calculate median values
heatmap_data_grouped <- heatmap_data_sorted %>%
  mutate(Age = age_sorted) %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)

# Extract final data for heatmap
age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped <- as.matrix(select(heatmap_data_grouped, -Age))

# Transpose data matrix (features as rows, age groups as columns)
heatmap_data_transposed <- t(heatmap_data_grouped)

# Separate cytokine and protein data
cytokine_colnames <- colnames(cytokine_common)
protein_colnames <- colnames(protein_common)
feature_types <- c(rep("Cytokine stimulus pair", length(cytokine_colnames)), 
                   rep("Protein", length(protein_colnames)))

cytokine_indices <- which(feature_types == "Cytokine stimulus pair")
protein_indices <- which(feature_types == "Protein")

cytokine_data <- heatmap_data_transposed[cytokine_indices, , drop = FALSE]
protein_data <- heatmap_data_transposed[protein_indices, , drop = FALSE]

# Z-score normalization function
zscore_data <- function(data_matrix) {
  zscore_matrix <- t(apply(data_matrix, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  return(zscore_matrix)
}

# Apply Z-score normalization to each feature type
cytokine_data_zscore <- zscore_data(cytokine_data)
protein_data_zscore <- zscore_data(protein_data)

# Age color mapping
age_col_fun <- colorRamp2(
  quantile(age_sorted, probs = seq(0, 1, length.out = 100)),
  colorRampPalette(c("#f2e8cf", "#a7c957"))(100)
)

# Calculate sample counts for each age group
age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))
age_counts_sorted <- age_counts[match(age_sorted, age_counts$Age), ]
sample_counts <- ifelse(is.na(age_counts_sorted$Count), 0, age_counts_sorted$Count)

# Create age labels with sample counts
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0, 
  #paste0(as.character(age_sorted), " (n=", sample_counts, ")"),
   as.character(age_sorted), 
  ""
)
# ========== Change 1: Add one line (Nature-style font, 6 pt) ==========
font_pt <- 5
row_height_mm <- 1.5 
heatmap_width_mm <- 45

##z-score
cytokine_zscore_range <- range(cytokine_data_zscore, na.rm = TRUE)
protein_zscore_range  <- range(protein_data_zscore,  na.rm = TRUE)

max_cytokine_abs <- max(abs(cytokine_zscore_range))
max_protein_abs  <- max(abs(protein_zscore_range))

cytokine_limit <- min(max_cytokine_abs, 2.5)
protein_limit  <- min(max_protein_abs,  2.5)

cytokine_col_fun <- colorRamp2(
  seq(-cytokine_limit, cytokine_limit, length.out = 100),
  colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)
protein_col_fun <- colorRamp2(
  seq(-protein_limit, protein_limit, length.out = 100),
  colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)


#age_labels[age_sorted == 70] <- "" #hide 70，overlap in the fig
# ========== Change 2: Add/update these in Heatmap() ==========
# Create cytokine heatmap
cytokine_heatmap <- Heatmap(
  cytokine_data_zscore,
  name = "Cytokine responds\nZ-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = font_pt),              # ADD
  column_names_gp = gpar(fontsize = font_pt),           # ADD
  col = cytokine_col_fun,
  row_dend_reorder = FALSE,
  row_dend_width = unit(5, "mm"),
    width = unit(heatmap_width_mm, "mm"),
    height = unit(nrow(cytokine_data_zscore) * row_height_mm, "mm"),
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun, height = unit(2, "mm")),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = font_pt), rot = 0),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = font_pt)
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#e63946", col = "white"),
      width = unit(2, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE,
    heatmap_legend_param = list(                          # ADD this block
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
      direction = "horizontal",
    legend_width = unit(15, "mm"),          
    grid_height = unit(2, "mm")    
)
)
# Create protein heatmap
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Proteomics\nZ-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt),              # ADD
  column_names_gp = gpar(fontsize = font_pt),           # ADD
  col = protein_col_fun,
  row_dend_reorder = FALSE,
    row_dend_width = unit(5, "mm"),
    width = unit(heatmap_width_mm, "mm"),
    height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#457b9d", col = "white"),
      width = unit(2, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE,
    heatmap_legend_param = list(                          # ADD this block
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
      direction = "horizontal",
    legend_width = unit(15, "mm"),          
    grid_height = unit(2, "mm")    
)
)
# ========== Change 3: Feature Type legend font ==========
feature_type_legend_separate <- Legend(
  title = "Feature Type",
  title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Cytokine\nresponds","Proteomics"),
  legend_gp = gpar(fill = c("#e63946","#457b9d")),
  labels_gp = gpar(fontsize = font_pt),
  grid_width = unit(2, "mm"),
    ncol = 2
)

# ========== Change 4: Sample count axis font ==========
sample_count_anno <- HeatmapAnnotation(
   `Sample size` = anno_barplot(
    sample_counts,
    gp = gpar(fill = "#a8dadc", col = NA, lwd = 0),
    height = unit(0.8, "cm"),
    axis_param = list(
      at = c(0, max(sample_counts)),
      labels = c("0", max(sample_counts)),
      gp = gpar(fontsize = font_pt)
    ),
    border = FALSE
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = font_pt)
)

heatmap_combined_separate <- sample_count_anno %v% (cytokine_heatmap %v% protein_heatmap)


# ========== Change 5: PDF size (Nature double-column) ==========
pdf("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_bar_zscore_spears_sig_feature_gender_nature.pdf", width = 3.62, height = 5.2, onefile = FALSE)

# draw() must be assigned so row_order() gets the initialized list (avoids warning)
ht_drawn <- draw(
  heatmap_combined_separate,
  column_title = "500FG",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  align_heatmap_legend = "heatmap_center",
  legend_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  legend_labels_gp = gpar(fontsize = font_pt),
  ht_gap = unit(c(2, 0.5), "mm")
)

# Export clustered row order from the drawn object (required for correct order)
row_orders <- row_order(ht_drawn)
cytokine_row_order <- rownames(cytokine_data_zscore)[row_orders[[1]]]
protein_row_order  <- rownames(protein_data_zscore)[row_orders[[2]]]
saveRDS(list(cytokine = cytokine_row_order, protein = protein_row_order),
        "/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_cyto_pro_row_gender_order.rds")
dev.off()

#29. **Display separate version on screen**
draw(
  heatmap_combined_separate,
  column_title = "500FG",                    
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
   annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  legend_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  legend_labels_gp = gpar(fontsize = font_pt),
  ht_gap = unit(c(2, 0.5), "mm")
)

save.image("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_bar_zscore_spears_sig_feature_gender_nature.Rdata")

###protein tissue link
protein_tissue <- read_excel("/vol/projects/yzhang/500FG_aging/input/olink/protein_cyto_sig_feature_tissue.xlsx")
protein_tissue <- protein_tissue[-c(1:20),]
head(protein_tissue)
dim(protein_tissue)
table(protein_tissue$GtexOrgan)

## Heatmap (tissue cluster) + Nature-style layout
# Prepare heatmap data (exclude Sample_ID and Age columns)
heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

# Sort data by age
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

# Group by age and calculate median values
heatmap_data_grouped <- heatmap_data_sorted %>%
  mutate(Age = age_sorted) %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)

# Extract final data for heatmap
age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped <- as.matrix(dplyr::select(heatmap_data_grouped, -Age))

# Transpose data matrix (features as rows, age groups as columns)
heatmap_data_transposed <- t(heatmap_data_grouped)

# Separate cytokine and protein data
cytokine_colnames <- colnames(cytokine_common)
protein_colnames <- colnames(protein_common)
feature_types <- c(
  rep("Cytokine stimulus pair", length(cytokine_colnames)),
  rep("Protein", length(protein_colnames))
)

cytokine_indices <- which(feature_types == "Cytokine stimulus pair")
protein_indices  <- which(feature_types == "Protein")

cytokine_data <- heatmap_data_transposed[cytokine_indices, , drop = FALSE]
protein_data  <- heatmap_data_transposed[protein_indices, , drop = FALSE]

# Z-score normalization
zscore_data <- function(data_matrix) {
  t(apply(data_matrix, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}

cytokine_data_zscore <- zscore_data(cytokine_data)
protein_data_zscore  <- zscore_data(protein_data)

# Age color mapping
age_col_fun <- circlize::colorRamp2(
  quantile(age_sorted, probs = seq(0, 1, length.out = 100), na.rm = TRUE),
  colorRampPalette(c("#f2e8cf", "#a7c957"))(100)
)

# Sample counts
age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))
age_counts_sorted <- age_counts[match(age_sorted, age_counts$Age), ]
sample_counts <- ifelse(is.na(age_counts_sorted$Count), 0, age_counts_sorted$Count)

# Age labels
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted, na.rm = TRUE), max(age_sorted, na.rm = TRUE)) | age_sorted %% 5 == 0,
  as.character(age_sorted),
  ""
)

## ---------- Nature-style size settings ----------
font_pt <- 5
row_height_mm <- 1.5
heatmap_width_mm <- 45

## ---------- Nature-style symmetric color scale ----------
cytokine_zscore_range <- range(cytokine_data_zscore, na.rm = TRUE)
protein_zscore_range  <- range(protein_data_zscore,  na.rm = TRUE)

cytokine_limit <- min(max(abs(cytokine_zscore_range)), 2.5)
protein_limit  <- min(max(abs(protein_zscore_range)),  2.5)

cytokine_col_fun <- circlize::colorRamp2(
  seq(-cytokine_limit, cytokine_limit, length.out = 100),
  colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)
protein_col_fun <- circlize::colorRamp2(
  seq(-protein_limit, protein_limit, length.out = 100),
  colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)

## ---------- Keep tissue cluster unchanged ----------
protein_rownames <- rownames(protein_data_zscore)
protein_organ_map <- protein_tissue$GtexOrgan
names(protein_organ_map) <- protein_tissue$Feature
protein_organs <- protein_organ_map[protein_rownames]
protein_organs[is.na(protein_organs) | protein_organs == "NA"] <- "Other"
          
# Cytokine heatmap
cytokine_heatmap <- ComplexHeatmap::Heatmap(
  cytokine_data_zscore,
  name = "Cytokine responds\nZ-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = font_pt),
  column_names_gp = grid::gpar(fontsize = font_pt),
  col = cytokine_col_fun,
  row_dend_reorder = FALSE,
  row_dend_width = grid::unit(5, "mm"),
  width = grid::unit(heatmap_width_mm, "mm"),
  height = grid::unit(nrow(cytokine_data_zscore) * row_height_mm, "mm"),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    Age = ComplexHeatmap::anno_simple(age_sorted, col = age_col_fun, height = grid::unit(2, "mm")),
    Age_Label = ComplexHeatmap::anno_text(age_labels, gp = grid::gpar(fontsize = font_pt), rot = 0),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = font_pt)
  ),
  left_annotation = ComplexHeatmap::rowAnnotation(
    Feature_Type = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "#e63946", col = "white"),
      width = grid::unit(2, "mm")
    )
  ),
  rect_gp = grid::gpar(col = NA),
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = font_pt),
    direction = "horizontal",
    legend_width = grid::unit(15, "mm"),
    grid_height = grid::unit(2, "mm")
  )
)

# Protein heatmap (tissue cluster unchanged)
protein_heatmap <- ComplexHeatmap::Heatmap(
  protein_data_zscore,
  name = "Proteomics\nZ-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = font_pt),
  column_names_gp = grid::gpar(fontsize = font_pt),
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  row_split = protein_organs,      # keep tissue clustering/splitting
  cluster_row_slices = FALSE,      # keep slice order stable
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = font_pt),
  row_dend_width = grid::unit(5, "mm"),
  row_gap = unit(0.2, "mm"), 
  width = grid::unit(heatmap_width_mm, "mm"),
  height = grid::unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  left_annotation = ComplexHeatmap::rowAnnotation(
    Feature_Type = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "#457b9d", col = "white"),
      width = grid::unit(2, "mm")
    )
  ),
  rect_gp = grid::gpar(col = NA),
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = font_pt),
    direction = "horizontal",
    legend_width = grid::unit(15, "mm"),
    grid_height = grid::unit(2, "mm")
  )
)

# Feature type legend
feature_type_legend_separate <- ComplexHeatmap::Legend(
  title = "Feature Type",
  title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Cytokine\nresponds", "Proteomics"),
  legend_gp = grid::gpar(fill = c("#e63946", "#457b9d")),
  labels_gp = grid::gpar(fontsize = font_pt),
  grid_width = grid::unit(2, "mm"),
  ncol = 2
)

# Sample count annotation
sample_count_anno <- ComplexHeatmap::HeatmapAnnotation(
  `Sample size` = ComplexHeatmap::anno_barplot(
    sample_counts,
    gp = grid::gpar(fill = "#a8dadc", col = NA, lwd = 0),
    height = grid::unit(0.8, "cm"),
    axis_param = list(
      at = c(0, max(sample_counts, na.rm = TRUE)),
      labels = c("0", max(sample_counts, na.rm = TRUE)),
      gp = grid::gpar(fontsize = font_pt)
    ),
    border = FALSE
  ),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(fontsize = font_pt)
)

# Combine heatmaps vertically
heatmap_combined_separate <- sample_count_anno %v% (cytokine_heatmap %v% protein_heatmap)

# Nature-size PDF + legends at bottom
pdf(
  "/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_spears_sig_feature_tissue_gender_nature.pdf",
  width = 4, height = 5.2, onefile = FALSE
)

ht_drawn <- ComplexHeatmap::draw(
  heatmap_combined_separate,
  column_title = "500FG",
  column_title_side = "top",
  column_title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  align_heatmap_legend = "heatmap_center",
  legend_title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
  legend_labels_gp = grid::gpar(fontsize = font_pt),
  ht_gap = grid::unit(c(2, 0), "mm")
)

dev.off()

# Display on screen
ComplexHeatmap::draw(
  heatmap_combined_separate,
  column_title = "500FG",
  column_title_side = "top",
  column_title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  legend_title_gp = grid::gpar(fontsize = font_pt, fontface = "bold"),
  legend_labels_gp = grid::gpar(fontsize = font_pt),
  ht_gap = grid::unit(c(2, 0.5), "mm")
)

save.image("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_bar_zscore_spears_sig_feature_tissue_gender.Rdata")






