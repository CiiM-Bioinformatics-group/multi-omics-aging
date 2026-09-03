library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)  # 用于自定义颜色
library(dendextend) 
library(NbClust)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(viridis)
library(readxl)

load("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_samefeature(500fg)_nature.Rdata")

load("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_samefeature(500fg)_nature_gender.Rdata")

# protein_500fg_names <- readRDS("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_protein_rownames.Rdata")
# print(protein_500fg_names)
# length(protein_500fg_names)
###set protein name same order 
ref <- readRDS("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_cyto_pro_row_gender_order.rds")
str(ref)
protein_500fg_names <- ref$protein
print(protein_500fg_names)
length(protein_500fg_names)

protein <- readRDS("/vol/projects/yzhang/300BCG/output/02_protein/olink_baseline_filter.RDS")
head(protein)
dim(protein)
sum(is.na(protein))

cytokine_corr <- read.csv("/vol/projects/yzhang/300BCG/output/01_cytokine/filtered_cytokine_name_replace_base.csv",row.names=1) ##no same feature，but use all cytokine feature
head(cytokine_corr)
dim(cytokine_corr)
cytokine_filtered <- cytokine_corr

#matching_columns <- intersect(colnames(protein), protein_500fg_names)
matching_columns <- protein_500fg_names[protein_500fg_names %in% colnames(protein)] ##change to same order
protein_subset <- protein[, matching_columns, drop = FALSE]
head(protein_subset)
dim(protein_subset)

pheno = read_excel("/vol/projects/CIIM/300BCG/201103 300BCG Age and sex data.xlsx")
head(pheno)
dim(pheno)
sum(is.na(pheno$Age))

# # 找到三个数据集共有的样本 ID
# common_samples <- intersect(rownames(protein_subset), pheno$PatientID)
# # 过滤出共同样本
# protein_filtered <- protein_subset[common_samples, , drop = FALSE]
# basicPhenos_filtered <- pheno[pheno$PatientID %in% common_samples, , drop = FALSE]
# basicPhenos_filtered <- as.data.frame(basicPhenos_filtered)
# # 重新设置 basicPhenos_filtered 的行名，使其与 common_samples 一致
# rownames(basicPhenos_filtered) <- basicPhenos_filtered$PatientID

# # 确保 basicPhenos_filtered 只包含共有样本
# basicPhenos_filtered <- basicPhenos_filtered[common_samples, , drop = FALSE]

# # 合并数据（保留所有特征）
# merged_data <- data.frame(
#   Sample_ID = common_samples,
#   Age = basicPhenos_filtered$Age,  # Age 对齐
#   protein_filtered  # 所有细胞因子特征
# )

# # 查看合并后的数据
# head(merged_data)
# dim(merged_data)  # 查看合并后的行数和列数


common_samples <- intersect(intersect(rownames(cytokine_filtered), 
                                      rownames(protein_filtered)), 
                            pheno$PatientID)

# Filter out common samples
cytokine_common <- cytokine_filtered[common_samples, , drop = FALSE]
protein_common <- protein_filtered[common_samples, , drop = FALSE]
pheno_common <- pheno[pheno$PatientID %in% common_samples, , drop = FALSE]

# Reset row names of basicPhenos_filtered to match common_samples
rownames(pheno_common) <- pheno_common$PatientID
# Ensure basicPhenos_filtered only contains common samples
pheno_common <- pheno_common[common_samples, , drop = FALSE]

# Merge data (retain all features)
merged_data <- data.frame(
  Sample_ID = common_samples,
  Age = pheno_common$Age,  # Age alignment
  cytokine_common,  # All cytokine features
  protein_common    # All protein features
)

# View merged data
head(merged_data)
dim(merged_data) 

write.csv(merged_data,"/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_merge_data.csv",row.names=TRUE)

heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

# 1. **Ensure Age is sorted in order**
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

# 2. **Remove NA and synchronize Age**
valid_rows <- complete.cases(heatmap_data_sorted)  # Ensure all rows have no NA
heatmap_data_sorted <- heatmap_data_sorted[valid_rows, , drop = FALSE]
age_sorted <- age_sorted[valid_rows]

# 3. **Convert DataFrame and merge Age**
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted)
heatmap_data_sorted$Age <- age_sorted

# 4. **Group by Age and calculate median**
heatmap_data_grouped <- heatmap_data_sorted %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)  # Ensure Age is still in ascending order

# Separate data BEFORE transposin
age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped_matrix <- as.matrix(select(heatmap_data_grouped, -Age))

# 7. **Create Feature Type vector** - All features are Protein type
# feature_types <- rep("Protein", ncol(heatmap_data_grouped))
# # 5. **Separate data by Feature Type BEFORE transposing**
# protein_cols <- which(feature_types == "Protein")

# protein_data_original <- heatmap_data_grouped_matrix[, protein_cols, drop = FALSE]

# # Transpose each data type separately
# protein_data <- t(protein_data_original)
protein_data_original <- heatmap_data_grouped_matrix
protein_data <- t(protein_data_original)

# 6. **Calculate sample counts for each Age group**
age_counts <- as.data.frame(table(merged_data$Age))  # Count original Age occurrences
colnames(age_counts) <- c("Age", "Count")  # Ensure clear column names
age_counts$Age <- as.numeric(as.character(age_counts$Age))  # Convert to numeric

# 7. **Merge sample counts to age_sorted (after grouping)**
age_counts_sorted <- age_counts[match(age_sorted, age_counts$Age), ]  # Ensure order matches
sample_counts <- ifelse(is.na(age_counts_sorted$Count), 0, age_counts_sorted$Count)  # Avoid NA

# 8. **Ensure sample_counts length matches heatmap columns**
# sample_counts should have the same length as the number of age groups in the heatmap
stopifnot(length(sample_counts) == length(age_sorted))

# 8. **Update Age Labels to show (n=sample_count)**
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0, 
  # paste0(as.character(age_sorted), " (n=", sample_counts, ")"),  # Format: "Age (n=X)" - commented out
  as.character(age_sorted),  # Show only age without sample count
  ""  # Keep blank elsewhere
)

# 23. **Apply z-score normalization to each Feature Type separately**
# Function to calculate z-score for each feature (row)
zscore_data <- function(data_matrix) {
  zscore_matrix <- t(apply(data_matrix, 1, function(x) {
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    
    # Handle constant values (zero variance)
    if (is.na(x_sd) || x_sd == 0) {
      return(rep(0, length(x)))
    }
    
    return((x - x_mean) / x_sd)
  }))
  return(zscore_matrix)
}

# Apply z-score normalization to each Feature Type
protein_data_zscore <- zscore_data(protein_data)

# Create appropriate color functions for z-score data
protein_zscore_range <- range(protein_data_zscore, na.rm = TRUE)

# Create symmetric color scales centered at 0
max_protein_abs <- max(abs(protein_zscore_range))

# Truncate to 2.5 standard deviations for better visualization
protein_limit <- min(max_protein_abs, 2.5)

# Color functions for z-score data with truncation
protein_col_fun <- colorRamp2(
  breaks = seq(-protein_limit, protein_limit, length.out = 100),
  colors = colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)

# Create breaks for legends (use truncated limits)
protein_breaks <- round(seq(-protein_limit, protein_limit, length.out = 7), 2)

age_col_fun <- colorRamp2(
  breaks = c(min(age_sorted), max(age_sorted)),
  colors = c("#f2e8cf", "#a7c957")
)

# Protein heatmap - 添加top_annotation
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Protein Expression",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = 9), rot = 90),
    annotation_name_side = "left"
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#457b9d", col = "white"),
      width = unit(4, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE
)

# 25. **Create legends for separate heatmaps**
protein_legend <- Legend(
  title = "Protein Expression", 
  col_fun = protein_col_fun,
  at = protein_breaks,
  labels = protein_breaks
)

feature_type_legend_separate <- Legend(
  title = "Feature Type",
  title_gp = gpar(fontsize = 12, fontface = "bold"),
#   at = c("Cytokine stimulus pair", "Protein"),
#   legend_gp = gpar(fill = c("#e63946", "#457b9d"))
# )
  at = c("Protein"),
  legend_gp = gpar(fill = c("#457b9d"))
)
# 26. **Create sample count bar chart annotation**
sample_count_anno <- HeatmapAnnotation(
  Sample_Count = anno_barplot(
    sample_counts,
    gp = gpar(fill = "#a8dadc", col = NA, lwd = 0),
    height = unit(3, "cm"),
    axis_param = list(
      at = c(0, max(sample_counts)),
      labels = c("0", max(sample_counts)),
      gp = gpar(fontsize = 8)
    ),
    border = FALSE  # Remove border
  ),
  annotation_name_side = "left"  # Move annotation name to left
)
# 27. **Combine separate heatmaps vertically with sample count annotation**
heatmap_combined_separate <- sample_count_anno %v% protein_heatmap

# 28. **Save separate version to PDF**
pdf("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_transposed_bar_zscore_cluster_samefeature(500fg).pdf", width = 18, height = 15, onefile = FALSE)

draw(
  heatmap_combined_separate,
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(15, 15, 15, 15), "mm")
)

dev.off()

#29. **Display separate version on screen**
draw(
  heatmap_combined_separate,
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)


save.image("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_samefeature(500fg).Rdata")

##same feature and same order
heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

# 1. **Ensure Age is sorted in order**
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

# 2. **Remove NA and synchronize Age**
valid_rows <- complete.cases(heatmap_data_sorted)  # Ensure all rows have no NA
heatmap_data_sorted <- heatmap_data_sorted[valid_rows, , drop = FALSE]
age_sorted <- age_sorted[valid_rows]

# 3. **Convert DataFrame and merge Age**
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted)
heatmap_data_sorted$Age <- age_sorted

# 4. **Group by Age and calculate median**
heatmap_data_grouped <- heatmap_data_sorted %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)  # Ensure Age is still in ascending order

# Separate data BEFORE transposin
age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped_matrix <- as.matrix(select(heatmap_data_grouped, -Age))

#7. **Create Feature Type vector** - All features are Protein type
feature_types <- rep("Protein", ncol(heatmap_data_grouped))
# 5. **Separate data by Feature Type BEFORE transposing**
protein_cols <- which(feature_types == "Protein")

protein_data_original <- heatmap_data_grouped_matrix[, protein_cols, drop = FALSE]

# # Transpose each data type separately
# protein_data <- t(protein_data_original)
protein_data_original <- heatmap_data_grouped_matrix
protein_data <- t(protein_data_original)

# 6. **Calculate sample counts for each Age group**
age_counts <- as.data.frame(table(merged_data$Age))  # Count original Age occurrences
colnames(age_counts) <- c("Age", "Count")  # Ensure clear column names
age_counts$Age <- as.numeric(as.character(age_counts$Age))  # Convert to numeric

# 7. **Merge sample counts to age_sorted (after grouping)**
age_counts_sorted <- age_counts[match(age_sorted, age_counts$Age), ]  # Ensure order matches
sample_counts <- ifelse(is.na(age_counts_sorted$Count), 0, age_counts_sorted$Count)  # Avoid NA

# 8. **Ensure sample_counts length matches heatmap columns**
# sample_counts should have the same length as the number of age groups in the heatmap
stopifnot(length(sample_counts) == length(age_sorted))

# 8. **Update Age Labels to show (n=sample_count)**
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0, 
  # paste0(as.character(age_sorted), " (n=", sample_counts, ")"),  # Format: "Age (n=X)" - commented out
  as.character(age_sorted),  # Show only age without sample count
  ""  # Keep blank elsewhere
)

# 23. **Apply z-score normalization to each Feature Type separately**
# Function to calculate z-score for each feature (row)
zscore_data <- function(data_matrix) {
  zscore_matrix <- t(apply(data_matrix, 1, function(x) {
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    
    # Handle constant values (zero variance)
    if (is.na(x_sd) || x_sd == 0) {
      return(rep(0, length(x)))
    }
    
    return((x - x_mean) / x_sd)
  }))
  return(zscore_matrix)
}

# Apply z-score normalization to each Feature Type
protein_data_zscore <- zscore_data(protein_data)

# Create appropriate color functions for z-score data
protein_zscore_range <- range(protein_data_zscore, na.rm = TRUE)

# Create symmetric color scales centered at 0
max_protein_abs <- max(abs(protein_zscore_range))

# Truncate to 2.5 standard deviations for better visualization
protein_limit <- min(max_protein_abs, 2.5)

# Color functions for z-score data with truncation
protein_col_fun <- colorRamp2(
  breaks = seq(-protein_limit, protein_limit, length.out = 100),
  colors = colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)

# Create breaks for legends (use truncated limits)
protein_breaks <- round(seq(-protein_limit, protein_limit, length.out = 7), 2)


age_col_fun <- colorRamp2(
  breaks = c(min(age_sorted), max(age_sorted)),
  colors = c("#f2e8cf", "#a7c957")
)

# Protein heatmap - 添加top_annotation
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Proteomics Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = 9), rot = 0),
    annotation_name_side = "left"
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#457b9d", col = "white"),
      width = unit(4, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE
)

# 25. **Create legends for separate heatmaps**
protein_legend <- Legend(
  title = "Proteomics Z-score", 
  col_fun = protein_col_fun,
  at = protein_breaks,
  labels = protein_breaks
)

feature_type_legend_separate <- Legend(
  title = "Feature Type",
  title_gp = gpar(fontsize = 12, fontface = "bold"),
#   at = c("Cytokine stimulus pair", "Protein"),
#   legend_gp = gpar(fill = c("#e63946", "#457b9d"))
# )
  at = c("Proteomics"),
  legend_gp = gpar(fill = c("#457b9d"))
)

# 26. **Create sample count bar chart annotation**
sample_count_anno <- HeatmapAnnotation(
  Sample_Count = anno_barplot(
    sample_counts,
    gp = gpar(fill = "#a8dadc", col = NA, lwd = 0),
    height = unit(3, "cm"),
    axis_param = list(
      at = c(0, max(sample_counts)),
      labels = c("0", max(sample_counts)),
      gp = gpar(fontsize = 8)
    ),
    border = FALSE  # Remove border
  ),
  annotation_name_side = "left"  # Move annotation name to left
)
# 27. **Combine separate heatmaps vertically with sample count annotation**
heatmap_combined_separate <- sample_count_anno %v% protein_heatmap


# 28. **Save separate version to PDF**
# pdf("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap__zscore_cluster_samefeature_order(500fg)_gender_cytokine.pdf", width = 18, height = 18, onefile = FALSE)

# draw(
#   heatmap_combined_separate,
#   show_heatmap_legend = TRUE,
#   show_annotation_legend = TRUE,
#    annotation_legend_list = list(feature_type_legend_separate),
#    heatmap_legend_side = "right",
#   annotation_legend_side = "right",
#   merge_legends = TRUE,
#   padding = unit(c(15, 15, 15, 15), "mm")
# )

# dev.off()

#29. **Display separate version on screen**
draw(
  heatmap_combined_separate,
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
  padding = unit(c(5, 5, 5, 5), "mm")
)


###nature style
# ========== Change 1: Add one line (Nature-style font, 6 pt) ==========
font_pt <- 5
row_height_mm <- 1.5 
heatmap_width_mm <- 45

age_labels[age_sorted == 70] <- "" #hide 70，overlap in the fig
# ========== Change 2: Add/update these in Heatmap() ==========
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Proteomics Z-score",
  height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt),              # ADD
  column_names_gp = gpar(fontsize = font_pt),           # ADD
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  width = unit(heatmap_width_mm, "mm"),
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun, height = unit(2, "mm")),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = font_pt), rot = 0),  
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = font_pt)
  ),
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
  at = c("Proteomics"),
  legend_gp = gpar(fill = c("#457b9d")),
  labels_gp = gpar(fontsize = font_pt),
    grid_height = unit(2, "mm"),    
  grid_width = unit(2, "mm")    
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

heatmap_combined_separate <- sample_count_anno %v% protein_heatmap

# ========== Change 5: PDF size (Nature double-column) ==========
pdf("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_zscore_cluster_samefeature_order(500fg)_nature_gender.pdf", width = 3.46, height = 3.5, onefile = FALSE)

draw(
  heatmap_combined_separate,
  column_title = "300BCG",                    
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
  legend_grid_height = unit(2, "mm"),
  legend_grid_width = unit(2, "mm"),
)

dev.off()

#29. **Display separate version on screen**
draw(
  heatmap_combined_separate,
  column_title = "300BCG",                    
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
   annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
    merge_legends = TRUE
)

###add cytokine version
# 0) basic matrix
heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

# 1) sort by age
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted_raw <- merged_data$Age[sorted_indices]

# 2) remove NA rows and sync age
valid_rows <- complete.cases(heatmap_data_sorted)
heatmap_data_sorted <- heatmap_data_sorted[valid_rows, , drop = FALSE]
age_sorted_raw <- age_sorted_raw[valid_rows]

# 3) group by age median (dplyr 1.1+ syntax)
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted)
heatmap_data_sorted$Age <- age_sorted_raw

heatmap_data_grouped <- heatmap_data_sorted %>%
  group_by(Age) %>%
  summarise(across(everything(), \(x) median(x, na.rm = TRUE))) %>%
  ungroup() %>%
  arrange(Age)

age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped_matrix <- as.matrix(dplyr::select(heatmap_data_grouped, -Age))

# 4) split cytokine/protein with fixed order
# cytokine order: keep current cytokine_common column order (your 7 cytokines)
cytokine_cols_use <- colnames(cytokine_common)[colnames(cytokine_common) %in% colnames(heatmap_data_grouped_matrix)]

# protein order: strictly follow 500fg reference order
protein_cols_use <- protein_500fg_names[protein_500fg_names %in% colnames(heatmap_data_grouped_matrix)]

cytokine_data_original <- heatmap_data_grouped_matrix[, cytokine_cols_use, drop = FALSE]
protein_data_original  <- heatmap_data_grouped_matrix[, protein_cols_use, drop = FALSE]

# transpose: rows=features, cols=age groups
cytokine_data <- t(cytokine_data_original)
protein_data  <- t(protein_data_original)

# 5) sample counts aligned to age_sorted
age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))

sample_counts <- age_counts$Count[match(age_sorted, age_counts$Age)]
sample_counts[is.na(sample_counts)] <- 0

# 6) age labels (recompute to avoid length mismatch)
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0,
  as.character(age_sorted),
  ""
)

# optional: hide one label if overlap
# age_labels[age_sorted == 70] <- ""

# 7) z-score function
zscore_data <- function(data_matrix) {
  zscore_matrix <- t(apply(data_matrix, 1, function(x) {
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    if (is.na(x_sd) || x_sd == 0) return(rep(0, length(x)))
    (x - x_mean) / x_sd
  }))
  return(zscore_matrix)
}

cytokine_data_zscore <- zscore_data(cytokine_data)
protein_data_zscore  <- zscore_data(protein_data)

# 8) color scales
cytokine_range <- range(cytokine_data_zscore, na.rm = TRUE)
protein_range  <- range(protein_data_zscore, na.rm = TRUE)

cytokine_limit <- min(max(abs(cytokine_range)), 2.5)
protein_limit  <- min(max(abs(protein_range)), 2.5)

cytokine_col_fun <- colorRamp2(
  breaks = seq(-cytokine_limit, cytokine_limit, length.out = 100),
  colors = colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)

protein_col_fun <- colorRamp2(
  breaks = seq(-protein_limit, protein_limit, length.out = 100),
  colors = colorRampPalette(c("#669bbc", "white", "#c1121f"))(100)
)

cytokine_breaks <- round(seq(-cytokine_limit, cytokine_limit, length.out = 7), 2)
protein_breaks  <- round(seq(-protein_limit,  protein_limit,  length.out = 7), 2)

age_col_fun <- colorRamp2(
  breaks = c(min(age_sorted), max(age_sorted)),
  colors = c("#f2e8cf", "#a7c957")
)

# 9) safety checks to prevent annotation length error
stopifnot(
  length(age_sorted) == ncol(cytokine_data_zscore),
  length(age_sorted) == ncol(protein_data_zscore),
  length(age_labels) == length(age_sorted),
  length(sample_counts) == length(age_sorted)
)

# 10) style
font_pt <- 8
row_height_mm <- 2.2
heatmap_width_mm <- 60

# 11) cytokine heatmap
cytokine_heatmap <- Heatmap(
  cytokine_data_zscore,
  name = "Cytokine Z-score",
  height = unit(nrow(cytokine_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = font_pt),
  col = cytokine_col_fun,
  row_dend_reorder = FALSE,
  width = unit(heatmap_width_mm, "mm"),
  top_annotation = HeatmapAnnotation(
    Age = anno_simple(age_sorted, col = age_col_fun, height = unit(2.5, "mm")),
    Age_Label = anno_text(age_labels, gp = gpar(fontsize = font_pt), rot = 0),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = font_pt)
  ),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#e63946", col = "white"),
      width = unit(3, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
    at = cytokine_breaks,
    direction = "horizontal",
    legend_width = unit(20, "mm"),
    grid_height = unit(2.5, "mm")
  )
)

# 12) protein heatmap
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Proteomics Z-score",
  height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt),
  column_names_gp = gpar(fontsize = font_pt),
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  width = unit(heatmap_width_mm, "mm"),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#457b9d", col = "white"),
      width = unit(3, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
    at = protein_breaks,
    direction = "horizontal",
    legend_width = unit(20, "mm"),
    grid_height = unit(2.5, "mm")
  )
)

# 13) legends
feature_type_legend_separate <- Legend(
  title = "Feature Type",
  title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Cytokine stimulus pair", "Proteomics"),
  legend_gp = gpar(fill = c("#e63946", "#457b9d")),
  labels_gp = gpar(fontsize = font_pt),
  grid_height = unit(2.5, "mm"),
  grid_width = unit(2.5, "mm")
)

# 14) sample count annotation
sample_count_anno <- HeatmapAnnotation(
  Sample_Count = anno_barplot(
    sample_counts,
    gp = gpar(fill = "#a8dadc", col = NA, lwd = 0),
    height = unit(1.0, "cm"),
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

# 15) combine
heatmap_combined_separate <- sample_count_anno %v% (cytokine_heatmap %v% protein_heatmap)

# 16) save PDF
pdf(
  "/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_zscore_cluster_samefeature_order(500fg)_nature_gender_cytokine.pdf",
  width = 6.5, height = 9.0, onefile = FALSE
)

draw(
  heatmap_combined_separate,
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(8, 8, 8, 8), "mm")
)

dev.off()

# 17) display in notebook
draw(
  heatmap_combined_separate,
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(4, 4, 4, 4), "mm")
)

### nature style ###add cytokine
font_pt <- 5
row_height_mm <- 1.5
heatmap_width_mm <- 45

age_labels[age_sorted == 70] <- ""

# cytokine heatmap
cytokine_heatmap <- Heatmap(
  cytokine_data_zscore,
  name = "Cytokine Expression",
  height = unit(nrow(cytokine_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = font_pt),
  col = cytokine_col_fun,
  row_dend_reorder = FALSE,
  width = unit(heatmap_width_mm, "mm"),
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
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
    direction = "horizontal",
    legend_width = unit(15, "mm"),
    grid_height = unit(2, "mm")
  )
)

# protein heatmap (unchanged logic)
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Proteomics Z-score",
  height = unit(nrow(protein_data_zscore) * row_height_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = font_pt),
  column_names_gp = gpar(fontsize = font_pt),
  col = protein_col_fun,
  row_dend_reorder = FALSE,
  width = unit(heatmap_width_mm, "mm"),
  left_annotation = rowAnnotation(
    Feature_Type = anno_block(
      gp = gpar(fill = "#457b9d", col = "white"),
      width = unit(2, "mm")
    )
  ),
  rect_gp = gpar(col = NA),
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_pt, fontface = "bold"),
    labels_gp = gpar(fontsize = font_pt),
    direction = "horizontal",
    legend_width = unit(15, "mm"),
    grid_height = unit(2, "mm")
  )
)

feature_type_legend_separate <- Legend(
  title = "Feature Type",
  title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  at = c("Cytokine responses", "Proteomics"),
  legend_gp = gpar(fill = c("#e63946", "#457b9d")),
  labels_gp = gpar(fontsize = font_pt),
  grid_height = unit(2, "mm"),
  grid_width = unit(2, "mm")
)

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

# combine: cytokine on top, protein below
heatmap_combined_separate <- sample_count_anno %v% (cytokine_heatmap %v% protein_heatmap)

pdf("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_zscore_samefeature_order(500fg)_nature_gender_cytokine.pdf",
    width = 3.46, height = 3.9, onefile = FALSE)

draw(
  heatmap_combined_separate,
  column_title = "300BCG",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)
dev.off()

#29. **Display separate version on screen**
draw(
  heatmap_combined_separate,
  column_title = "300BCG",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
  annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)

save.image("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_samefeature(500fg)_nature_gender.Rdata")

save.image("/vol/projects/yzhang/300BCG/output/02_protein/linear_age/300bcg_heatmap_samefeature(500fg)_nature_gender_cytokine.Rdata")

out_300bcg <- list(
  protein_data = protein_data,
  sample_counts       = sample_counts,
  age_sorted          = age_sorted,
  age_labels          = age_labels,
  age_col_fun         = age_col_fun
)
saveRDS(out_300bcg, "/vol/projects/yzhang/300BCG/output/02_protein/linear_age/heatmap_300bcg_nature.rds")
