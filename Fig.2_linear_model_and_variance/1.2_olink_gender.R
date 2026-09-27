# Olink (proteomics) ~ Age association analysis (adjusted for Gender)

library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(networkD3)
library(tibble)
library(ggalluvial)
library(RColorBrewer)
library(scales)
library(viridis)

source("utils.R")

layer_output_dir <- file.path(project_dir, "output", "02_Olink")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

olink <- readRDS(file.path(raw_data_dir, "olinkData", "NPX_FG500.RDS"))

basicPhenos <- read.csv(basic_phenos_path, row.names = 1)
rownames(basicPhenos) <- basicPhenos$ID_500fg

idx <- intersect(rownames(olink), basicPhenos$ID_500fg)  # 458 shared samples
olink_filter <- olink[rownames(olink) %in% idx, , drop = FALSE]
basicPhenos_filter_match <- basicPhenos %>%
  filter(ID_500fg %in% rownames(olink_filter)) %>%
  slice(match(rownames(olink_filter), ID_500fg))
basicPhenos_filter_match <- clean_gender(basicPhenos_filter_match)

# model 1: age + protein (age not scaled) + Gender
final_result <- run_age_gender_lm(
  olink_filter,
  basicPhenos_filter_match,
  transform = function(x) log(x + 1)
)

write.csv(final_result, file = file.path(layer_output_dir, "olink_age_liner_gender_FDR.csv"))

# 显著性分类（上调/下调/非显著），用于火山图着色
final_result <- final_result %>%
  mutate(
    sig = ifelse(p < 0.05 & estimate > 0, "upregulated",
                 ifelse(p < 0.05 & estimate < 0, "downregulated", "not significant")),
    log_p = -log10(p)
  )

top4_combined <- bind_rows(
  final_result %>% filter(sig == "upregulated") %>% arrange(desc(log_p)) %>% slice(1:3),
  final_result %>% filter(sig == "downregulated") %>% arrange(desc(log_p)) %>% slice(1:3)
)

ppi <- 300
png(file.path(layer_output_dir, "volcano_age_olink_liner_top3(2).png"),
    width = 5 * ppi, height = 4 * ppi, res = ppi)

ggplot(final_result, aes(x = estimate, y = log_p)) +
  geom_point(aes(color = sig), size = 3, alpha = 0.6) +
  scale_color_manual(values = c(
    "upregulated" = "#FFB3E6",
    "downregulated" = "#80B1D3",
    "not significant" = "grey"
  )) +
  theme_classic() +
  xlab("Estimate (Effect Size)") +
  ylab("-log10(p)") +
  ggtitle("Proteomics ~ Age Volcano Plot") +
  annotate(
    "text",
    x = max(final_result$estimate, na.rm = TRUE) * -1.3,
    y = max(final_result$log_p, na.rm = TRUE) - 3,
    label = paste0(
      "Upregulated: ", sum(final_result$sig == "upregulated"),
      " (", round(mean(final_result$sig == "upregulated") * 100, 2), "%)",
      "\nDownregulated: ", sum(final_result$sig == "downregulated"),
      " (", round(mean(final_result$sig == "downregulated") * 100, 2), "%)"
    ),
    size = 4, hjust = 0, color = "black"
  ) +
  geom_text_repel(
    data = top4_combined,
    aes(x = estimate, y = log_p, label = olink),
    size = 3.5, color = "black", box.padding = 0.5, point.padding = 0.5
  )

dev.off()

write_sig_results(final_result, file.path(layer_output_dir, "olink_age_liner_gender_FDR_sig.csv"))
