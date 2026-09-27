# Cytokine ~ Age association analysis (adjusted for Gender)
# Run from the repository root. No data is included in this repo -- see
# README for the expected input layout.

library(readxl)
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
library(readr)

source("utils.R")

layer_output_dir <- file.path(project_dir, "output", "01_cytokine")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

cytokine_path <- file.path(project_dir, "input", "cytokine", "filtered_cytokine_filled_trans.csv")
cytokine <- read.csv(cytokine_path)
rownames(cytokine) <- cytokine[, 1]
cytokine <- cytokine[, -1]  # 489 samples x 91 cytokines

basicPhenos <- read.csv(basic_phenos_path, row.names = 1)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(cytokine, basicPhenos)
cytokine_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 1: age + cytokine (age not scaled) + Gender
final_result <- run_age_gender_lm(
  cytokine_filter,
  basicPhenos_filter_match,
  transform = function(x) log(x + 1)
)

write.csv(final_result, file = file.path(layer_output_dir, "cytokine_age_liner_gender_FDR.csv"))

write_sig_results(
  final_result,
  file.path(layer_output_dir, "cytokine_age_liner_gender_FDR_sig.csv")
)
