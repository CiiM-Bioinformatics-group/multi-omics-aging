library(dplyr)
library(magrittr)
library(caret)
library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(magrittr)
library(reshape2)
library(caret)
library(tidyr)
library(ggrepel)
library(scales)
library(tidyverse)
library(readr)
library(glue)
library(jsonlite)
library(parallel)

library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)  # 用于自定义颜色
library(dendextend) 
library(NbClust)
library(ggpubr)
library(viridis)

load("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/heatmap_prime_zscore_same_feature(500fg).Rdata")

#olink <- readRDS("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/input/prime_olink_common.rds")
#cytokine <- readRDS("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/input/prime_cytokine_common.rds")
phenotype <- readRDS("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/input/prime_phenotype_common.rds")

#head(cytokine)
#dim(cytokine)
#head(olink)
#dim(olink)
head(phenotype)
dim(phenotype)

###raw data 
olink <- readRDS("/vol/projects/yzhang/BCG_prime/03_protein/input/olink_raw_cleaned.RDS")
head(olink)
dim(olink)
olink_common <- olink[rownames(olink) %in% rownames(phenotype), ]
olink_common <- olink_common[match(rownames(phenotype), rownames(olink_common)), ] ###change order 
head(olink_common)
dim(olink_common)

cytokine <- readRDS("/vol/projects/yzhang/BCG_prime/02_cytokine/input/cytokine_trans_raw.RDS")
head(cytokine)
dim(cytokine)
cytokine_common <- cytokine[rownames(cytokine) %in% rownames(phenotype), ]
cytokine_common <- cytokine_common[match(rownames(phenotype), rownames(cytokine_common)), ]
head(cytokine_common)
dim(cytokine_common)

cytokine_common_filled <- apply(cytokine_common, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
olink_common_filled <- apply(olink_common, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))    
head(cytokine_common_filled)
head(olink_common_filled)

# cytokine_500fg_rownames <- readRDS("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_cytokine_rownames.Rdata")
# protein_500fg_rownames <- readRDS("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_protein_rownames.Rdata")
ref <- readRDS("/vol/projects/yzhang/500FG_aging/output/14_AA/heatmap_500fg_cyto_pro_row_gender_order.rds")
str(ref)
protein_500fg_names <- ref$protein
print(protein_500fg_names)
length(protein_500fg_names)

#matching_columns <- intersect(colnames(olink_common_filled), protein_500fg_rownames)
matching_columns <- protein_500fg_names[protein_500fg_names %in% colnames(olink_common_filled)]
protein_subset    <- olink_common_filled[, matching_columns, drop = FALSE]
head(protein_subset)
dim(protein_subset)

cytokine_mapping <- read.csv("/vol/projects/yzhang/BCG_prime/02_cytokine/input/cytokine_prime_name_mapping.csv")
head(cytokine_mapping)
colnames(cytokine_common_filled) <- cytokine_mapping$mapping_name[match(colnames(cytokine_common_filled), cytokine_mapping$original_name)]
head(cytokine_common_filled)

colnames(cytokine_common_filled)
cytokine_500fg_rownames

matching_columns <- intersect(colnames(cytokine_common_filled), cytokine_500fg_rownames)
cytokine_subset <- cytokine_common_filled[, matching_columns, drop = FALSE]
head(cytokine_subset)
dim(cytokine_subset)

# 合并数据（保留所有特征）
merged_data <- data.frame(
  Sample_ID = rownames(phenotype),
  Age = phenotype$Age,  # Age 对齐
  protein_subset    # 所有蛋白质特征
)

# 查看合并后的数据
head(merged_data)
dim(merged_data) 


heatmap_data <- merged_data[, !(colnames(merged_data) %in% c("Sample_ID", "Age"))]

write.csv(merged_data,"/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/BCG_prime_heatmap_merged_data.csv",row.names=TRUE)

# 1. **Ensure Age is sorted in order**
sorted_indices <- order(merged_data$Age, na.last = FALSE)
heatmap_data_sorted <- heatmap_data[sorted_indices, , drop = FALSE]
age_sorted <- merged_data$Age[sorted_indices]

# 2. **Convert DataFrame and merge Age**
heatmap_data_sorted <- as.data.frame(heatmap_data_sorted)
heatmap_data_sorted$Age <- age_sorted

# 3. **Group by Age and calculate median**
heatmap_data_grouped <- heatmap_data_sorted %>%
  group_by(Age) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Age)

# 4. **Prepare age annotation color**
age_col_fun <- colorRamp2(
  quantile(age_sorted, probs = seq(0, 1, length.out = 100)),
  colorRampPalette(c("#f2e8cf", "#a7c957"))(100)
)

# 5. **Extract final data**
age_sorted <- heatmap_data_grouped$Age
heatmap_data_grouped <- as.matrix(select(heatmap_data_grouped, -Age))

# 6. **Calculate sample counts for each Age group**
age_counts <- as.data.frame(table(merged_data$Age))
colnames(age_counts) <- c("Age", "Count")
age_counts$Age <- as.numeric(as.character(age_counts$Age))

# 7. **Merge sample counts to age_sorted**
age_counts_sorted <- age_counts[match(age_sorted, age_counts$Age), ]
sample_counts <- ifelse(is.na(age_counts_sorted$Count), 0, age_counts_sorted$Count)

# 8. **Update Age Labels to show (n=sample_count)**
age_labels <- ifelse(
  age_sorted %in% c(min(age_sorted), max(age_sorted)) | age_sorted %% 5 == 0, 
  paste0(as.character(age_sorted), " (n=", sample_counts, ")"),
  ""
)

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
  colors = c("#cfd796", "#a7c957")
)

# Protein heatmap - 添加top_annotation
protein_heatmap <- Heatmap(
  protein_data_zscore,
  name = "Protein Expression",
  #cluster_rows = TRUE, change to same order
    cluster_rows = FALSE,
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

#28. **Save separate version to PDF**
pdf("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/heatmap_prime_zscore_same_order_feature(500fg)_gender.pdf", width = 18, height = 15, onefile = FALSE)

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


###nature style
# ========== Change 1: Add one line (Nature-style font, 6 pt) ==========
font_pt <- 5
row_height_mm <- 1.5 
heatmap_width_mm <- 45

#age_labels[age_sorted == 70] <- "" #hide 70，overlap in the fig
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
pdf("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/heatmap_prime_zscore_same_feature(500fg)_nature_gender.pdf", width = 3.46, height = 3.5, onefile = FALSE)

draw(
  heatmap_combined_separate,
  column_title = "BCG-Prime",                    
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
  column_title = "BCG-Prime",                    
  column_title_side = "top",
  column_title_gp = gpar(fontsize = font_pt, fontface = "bold"),
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE,
   annotation_legend_list = list(feature_type_legend_separate),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
    merge_legends = TRUE
)

save.image("/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/heatmap_prime_zscore_same_feature(500fg)_gender.Rdata")

out_prime <- list(
  protein_data = protein_data,
  sample_counts      = sample_counts,
  age_sorted         = age_sorted,
  age_labels         = age_labels,
  age_col_fun        = age_col_fun
)
saveRDS(out_prime, "/vol/projects/yzhang/BCG_prime/05_multi-omics_model/output/heatmap/heatmap_prime_nature.rds")
