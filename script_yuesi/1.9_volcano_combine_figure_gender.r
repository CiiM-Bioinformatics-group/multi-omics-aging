library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)

 #age not scale; correct gender
olink <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv",row.names=1) 
cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv",row.names=1)
metabolite <- read.csv("/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR.csv",row.names=1) 
methylation <- readRDS("/vol/projects/yzhang/500FG_aging/output/05_methylation/methylation_age_liner_gender_FDR.RDS") 
# rownames(methylation)<-methylation$X
# methylation$X <- NULL
methylation <- methylation %>% rename(feature = Mvalue)
hormone <- read.csv("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR.csv",row.names=1)
cellcounts <- read.csv("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv",row.names=1)  
immunoglobulin <- read.csv("/vol/projects/yzhang/500FG_aging/output/08_immunoglobulin/immunoglobulin_age_gender_FDR.csv",row.names=1)
microbiome <- read.csv("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR.csv",row.names=1)
platelet<- read.csv("/vol/projects/yzhang/500FG_aging/output/09_platelet/platelet_age_gender_FDR.csv",row.names=1)

head(olink)
head(cytokine)
head(metabolite)
head(methylation)
head(cellcounts)
head(microbiome)
head(hormone)
head(immunoglobulin)
head(platelet)

##nature style ###small size
# ------------------------------------------------------------------------------
# Add significance classification and -log10(padj) column
# ------------------------------------------------------------------------------
process_data <- function(data, level_name) {
  data %>%
    mutate(
      sig = ifelse(padj < 0.05, "significant", "not significant"),
      sig_direction = case_when(
        padj < 0.05 & estimate > 0 ~ "positive correlation",
        padj < 0.05 & estimate < 0 ~ "negative correlation",
        TRUE ~ "not significant"
      ),
      log_padj = -log10(padj),
      level = level_name
    )
}

# Layer display names (used for legend and labels)
titles <- c("DNA methylation", "Proteomics", "Metabolomics", "Cytokine response",
            "Immune cell counts", "Microbiomes", "Circulating endocrine traits", "Circulating immune traits")
# Datasets per layer: order must match titles 1:1
datasets <- list(
  methylation,    # DNA methylation
  olink,          # Proteomics
  metabolite,     # Metabolomics
  cytokine,       # Cytokine response
  cellcounts,     # Immune cell counts
  microbiome,     # Microbiomes
  hormone,        # Circulating endocrine traits
  immunoglobulin  # Circulating immune traits
)
names(datasets) <- titles

# Process and combine all layers
final_result_combined <- bind_rows(
  lapply(names(datasets), function(level_name) process_data(datasets[[level_name]], level_name))
)

# Group variable for color (significant vs not significant)
final_result_combined <- final_result_combined %>%
  mutate(level_sig = ifelse(sig == "not significant", "not significant", level))
color_palette <- setNames(
  # c("#B2DF8A", "#6A3D9A", "#FFFF99", "#FF7F00", "#A6CEE3", "#CAB2D6", "#e8d5b7", "#dcd6f7"),
     c("#B3DE69", "#FB8072", "#FFFFB3", "#80B1D3", "#BEBADA", "#FCCDE5", "#FDB462", "#8DD3C7"), 
  titles
)
legend_order <- c(titles, "not significant")
legend_colors <- c(color_palette, "not significant" = "grey")
colors <- final_result_combined %>%
  distinct(level_sig) %>%
  mutate(
    color = if_else(level_sig == "not significant", "grey", unname(color_palette[level_sig]))
  )
# Put Proteomics and Cytokine response on top layer: reorder data so they are drawn last
plot_data <- final_result_combined %>%
  mutate(.layer = case_when(
    level_sig %in% c("Proteomics", "Cytokine response") ~ 2L,
    TRUE ~ 1L
  )) %>%
  arrange(.layer)

# To add a title, uncomment: + ggtitle("Your title")
p <- ggplot(plot_data, aes(x = estimate, y = log_padj)) +
  #geom_point(aes(color = level_sig), size = 2, alpha = 0.9) +
  geom_point_rast(aes(color = level_sig), size = 2, alpha = 0.9, raster.dpi = 600)+ ###change the format of dots drawing
  scale_color_manual(values = legend_colors, breaks = legend_order) +
  scale_x_continuous(
    "Estimate (Effect Size)",
    expand = expansion(mult = 0.05)
  ) +
  scale_y_continuous(
    expression(-log[10](P[adj])),
    limits = c(-4, NA),
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::extended_breaks(n = 8)
  ) +
theme_classic(base_size = 5) +
  theme(
    axis.title = element_text(size = 5, colour = "black"),
    axis.text = element_text(size = 5, colour = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks = element_line(size = 0.5, colour = "black"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 5),
    legend.spacing.y = unit(2, "pt"),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

# 2) Replace the two ggsave + print lines with:
fig_w_mm <- 88
fig_h_mm <- 88 * 7 / 10  # same aspect as 10x7 in
ggsave("/vol/projects/yzhang/500FG_aging/output/13_volcano_combined/multi_level_volcano_gender_complete_nature_small.pdf", plot = p, width = fig_w_mm, height = fig_h_mm, units = "mm")
ggsave("/vol/projects/yzhang/500FG_aging/output/13_volcano_combined/multi_level_volcano_gender_complete_nature_small.png", plot = p, width = fig_w_mm, height = fig_h_mm, units = "mm", dpi = 300)
print(p)


save.image("/vol/projects/yzhang/500FG_aging/output/13_volcano_combined/multi_level_volcano_gender_complete.Rdata")
