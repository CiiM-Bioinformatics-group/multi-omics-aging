library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(patchwork)
library(ggpubr)

load("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/age_trajectory_nature.Rdata")

cellcount_cor <- read.csv("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv",row.names=1)
head(cellcount_cor)
dim(cellcount_cor)
cellcount_sig <- cellcount_cor[cellcount_cor$sig=="sig",]
head(cellcount_sig)
dim(cellcount_sig)

cellcounts <- read.csv("/vol/projects/yzhang/500FG_aging/input/cellCounts/cellcounts_name_replace.csv",row.names=1,check.names = FALSE)
head(cellcounts)
dim(cellcounts)
cellcounts_info <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_cellcounts_info.txt")
head(cellcounts_info)
colnames(cellcounts) <- cellcounts_info$finalName[match(colnames(cellcounts), cellcounts_info$trait)]
head(cellcounts)
dim(cellcounts)

saveRDS(cellcounts,"/vol/projects/yzhang/500FG_aging/input/cellCounts/cellcounts_finalname_replace.rds")

cellcount_exp <- cellcounts[, colnames(cellcounts) %in% cellcount_sig$feature]
head(cellcount_exp)
dim(cellcount_exp)

pheno <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv",row.names=1)
rownames(pheno) <- pheno$ID_500fg
head(pheno)

mat <- as.matrix(t(cellcount_exp))
# Scale by row (each cell type scaled across individuals)
mat_scaled <- t(scale(t(mat)))

# ------------------------------------------------------------------------------
# 2. Cluster columns (individuals) -> Group 1, 2, 3 (immunotype)
# ------------------------------------------------------------------------------
hc_col <- hclust(dist(t(mat_scaled), method = "euclidean"), method = "ward.D2")
immunotype <- cutree(hc_col, k = 3)
group_label <- factor(paste0("Group ", immunotype[colnames(mat)]), levels = paste0("Group ", 1:3))

# ------------------------------------------------------------------------------
# 3. Cluster rows (cell types) -> C1, C2, C3
# ------------------------------------------------------------------------------
hc_row <- hclust(dist(mat_scaled, method = "euclidean"), method = "ward.D2")
feature_cluster <- cutree(hc_row, k = 3)

# ------------------------------------------------------------------------------
# 4. Heatmap colors: scaled proportion -5 to 5 (blue-white-red)
# ------------------------------------------------------------------------------
col_fun <- colorRamp2(
  seq(-5, 5, length.out = 9),
  rev(c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"))
)

# ------------------------------------------------------------------------------
# 5. Draw heatmap (rows = features, cols = individuals); Group 1/2/3 as column titles
# ------------------------------------------------------------------------------
ht <- Heatmap(
  mat_scaled,
  name = "Scaled cell proportion",
  cluster_rows = hc_row,
  cluster_columns = hc_col,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = 3,
  column_split = 3,
  column_title = c("Group 1", "Group 2", "Group 3"),
  column_title_gp = gpar(fontsize = 5, fontface = "bold"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  # show_row_names = FALSE,
  row_title = c("Cluster 1", "Cluster 2", "Cluster 3"),
  row_title_gp = gpar(fontsize = 5, fontface = "bold"),
  col = col_fun,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    direction = "horizontal",
    legend_width = unit(1.5, "cm"),
    grid_height = unit(1, "mm") 
  )
)

# 62 mm -> inches (1 inch = 25.4 mm)
pdf("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/cellcount_age_trajectory_heatmap_gender.pdf", width = 88/25.4, height = 80/25.4)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

draw(ht, newpage = TRUE)

top_anno <- HeatmapAnnotation(
  Age = anno_simple(
    age_vec,
    col = colorRamp2(
      range(age_vec, na.rm = TRUE),
      c("#f7f7f7", "#2166ac")
    )
  ),
  annotation_name_gp = gpar(fontsize = 5)
  # no annotation_legend_param
)

# Add top_annotation in Heatmap()
ht <- Heatmap(
  mat_scaled,
  name = "Scaled cell proportion",
  top_annotation = top_anno,   # add this line
  cluster_rows = hc_row,
  cluster_columns = hc_col,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = 3,
  column_split = 3,
  column_title = c("Group 1", "Group 2", "Group 3"),
  column_title_gp = gpar(fontsize = 5, fontface = "bold"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  # show_row_names = FALSE,
  row_title = c("Cluster 1", "Cluster 2", "Cluster 3"),
  row_title_gp = gpar(fontsize = 5, fontface = "bold"),
  col = col_fun,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    direction = "horizontal",
    legend_width = unit(1.5, "cm"),
    grid_height = unit(1, "mm") 
  )
)

pdf("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/cellcount_age_trajectory_heatmap_age.pdf", width = 88/25.4, height = 80/25.4)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

draw(ht, newpage = TRUE)

df_age <- data.frame(
  age   = age_vec,
  group = group_label
)

aggregate(age ~ group, data = df_age, FUN = function(x) c(
  n = length(x),
  mean = mean(x, na.rm = TRUE),
  sd   = sd(x, na.rm = TRUE),
  median = median(x, na.rm = TRUE),
  IQR   = IQR(x, na.rm = TRUE)
))


# One-way ANOVA (assumes normality and similar variances)
fit_aov <- aov(age ~ group, data = df_age)
summary(fit_aov)
kruskal.test(age ~ group, data = df_age)

# 3. Pairwise comparisons (optional)
# Wilcoxon rank-sum test for each pair, with BH correction
pairwise.wilcox.test(df_age$age, df_age$group, p.adjust.method = "BH", exact = FALSE)


p <- ggplot(df_age, aes(x = group, y = age, fill = group)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.6, alpha = 0.9) +
  scale_fill_manual(values = c("Group 1" = "#92c5de", "Group 2" = "#f4a582", "Group 3" = "#d6604d")) +
  labs(
    x = "Group",
    y = "Age (years)"
  ) +
  theme_classic(base_size = 5) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(color = "black", size = 5),
    axis.ticks = element_line(linewidth = 0.25)
  )

p_sig <- p + stat_compare_means(
  comparisons = list(
    c("Group 1", "Group 2"),
    c("Group 1", "Group 3"),
    c("Group 2", "Group 3")
  ),
  method = "wilcox.test",
  label = "p.signif",
  tip.length = 0.02,
  size = 2.2   # bracket/label size relative to base (adjust if asterisks too small)
)

# Save PDF: 80 mm x 42 mm (1 inch = 25.4 mm)
ggsave(
  "/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/age_by_group_violin2_gender.pdf",
  plot = p_sig,
  width = 70 / 25.4,
  height = 47 / 25.4,
  units = "in",
  device = "pdf"
)
print(p)

# Feature cluster: map cell type -> C1/C2/C3 (from heatmap row clustering; rownames(mat) = cell types)
fc_vec_feature <- setNames(paste0("C", feature_cluster), rownames(mat))

# Raw proportion matrix: rows = individuals, cols = cell types (use cellcount_exp, not mat)
mat_raw <- as.matrix(cellcount_exp)

# Per-individual mean proportion per feature cluster.
# If cellcount_exp is in 0-1 (proportion), use *100 to get %; if already 0-100 (%), do not multiply.
# Check: if max value in data > 1.5, assume already in %.
already_pct <- max(mat_raw, na.rm = TRUE) > 1.5
cluster_means <- mat_raw %>%
  as.data.frame() %>%
  tibble::rownames_to_column("individual") %>%
  pivot_longer(-individual, names_to = "cell_type", values_to = "proportion") %>%
  mutate(feature_cluster = fc_vec_feature[.data$cell_type]) %>%
  filter(!is.na(.data$feature_cluster)) %>%
  group_by(individual, feature_cluster) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE) * if (already_pct) 1 else 100, .groups = "drop")

# Add age; drop missing age
age_vec_traj <- setNames(pheno[rownames(mat_raw), "Age"], rownames(mat_raw))
traj_df <- cluster_means %>%
  mutate(Age = age_vec_traj[.data$individual]) %>%
  filter(is.finite(.data$Age)) %>%
  mutate(feature_cluster = factor(.data$feature_cluster, levels = c("C1", "C2", "C3")))

# Nature-style colors (colorblind-friendly, distinct)
# cluster_colors <- c("C1" = "#0173b2", "C2" = "#029e73", "C3" = "#de8f05")
 # cluster_colors <- c("C1" = "#0EA5E9", "C2" = "#10B981", "C3" = "#F59E0B")
# cluster_colors <- c("C1" = "#0072B2", "C2" = "#009E73", "C3" = "#D55E00")
cluster_colors <- c("C1" = "#8DD3C7", "C2" = "#80B1D3", "C3" = "#FB8072")
# One panel: 3 LOESS curves (C1, C2, C3 vs age)
# Y-axis: "Mean cell proportion (%)" = mean proportion of cells in that cluster, in percent (correct)
p_traj <- ggplot(traj_df, aes(x = Age, y = mean_prop, color = feature_cluster, fill = feature_cluster)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = cluster_colors, name = "Feature cluster") +
  scale_fill_manual(values = cluster_colors, name = "Feature cluster", guide = "none") +
  labs(x = "Age (years)", y = "Mean cell proportion (%)") +
  theme_minimal(base_size = 5) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )
print(p_traj)            
# 85 mm -> inches
ggsave("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/age_trajectory_C1_C2_C3.pdf", p_traj, width = 85/25.4, height = 60/25.4)


# Feature cluster: map cell type -> C1/C2/C3 (from heatmap row clustering; rownames(mat) = cell types)
fc_vec_feature <- setNames(paste0("C", feature_cluster), rownames(mat))

# Raw proportion matrix: rows = individuals, cols = cell types (use cellcount_exp, not mat)
mat_raw <- as.matrix(cellcount_exp)

# Per-individual mean proportion per feature cluster.
# If cellcount_exp is in 0-1 (proportion), use *100 to get %; if already 0-100 (%), do not multiply.
# Check: if max value in data > 1.5, assume already in %.
already_pct <- max(mat_raw, na.rm = TRUE) > 1.5
cluster_means <- mat_raw %>%
  as.data.frame() %>%
  tibble::rownames_to_column("individual") %>%
  pivot_longer(-individual, names_to = "cell_type", values_to = "proportion") %>%
  mutate(feature_cluster = fc_vec_feature[.data$cell_type]) %>%
  filter(!is.na(.data$feature_cluster)) %>%
  group_by(individual, feature_cluster) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE) * if (already_pct) 1 else 100, .groups = "drop")
# Nature-style colors (colorblind-friendly, distinct)
cluster_colors <- c("C1" = "#8DD3C7", "C2" = "#80B1D3", "C3" = "#FB8072")
# Add age; drop missing age
age_vec_traj <- setNames(pheno[rownames(mat_raw), "Age"], rownames(mat_raw))
traj_df <- cluster_means %>%
  mutate(Age = age_vec_traj[.data$individual]) %>%
  filter(is.finite(.data$Age)) %>%
  mutate(feature_cluster = factor(.data$feature_cluster, levels = c("C1", "C2", "C3")))
            
p_traj <- ggplot(traj_df, aes(x = Age, y = mean_prop, color = feature_cluster, fill = feature_cluster)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = cluster_colors, name = "Feature cluster") +
  scale_fill_manual(values = cluster_colors, name = "Feature cluster", guide = "none") +
  labs(x = "Age (years)", y = "Mean cell proportion (%)") +
  theme_minimal(base_size = 5) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )
print(p_traj)
ggsave("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/cellcount_age_trajectory_gender.pdf", p_traj, width = 85/25.4, height = 42/25.4)

save.image("/vol/projects/yzhang/500FG_aging/output/20_cellcount_age_trajectory/age_trajectory_nature_gender.Rdata")
