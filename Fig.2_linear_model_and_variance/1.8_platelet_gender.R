# Platelet count ~ Age association analysis (adjusted for Gender)
# Run from the repository root. No data is included in this repo -- see
# README for the expected input layout.
#
# NOTE (see README "Known issues"): unlike the other layers, this script
# stops after writing the main FDR CSV (no RData snapshot, no sig-only
# export), matching the original script exactly.

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

layer_output_dir <- file.path(project_dir, "output", "09_platelet")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

platelet <- read.table(
  file.path(raw_data_dir, "Molecular_phenotype", "500FG_log2_platelet_count.txt"),
  header = TRUE, sep = "", stringsAsFactors = FALSE
)

basicPhenos <- read.csv(basic_phenos_path)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(platelet, basicPhenos)
platelet_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 2: age + platelet using FDR, age not scaled + gender
final_result <- run_age_gender_lm(
  platelet_filter,
  basicPhenos_filter_match,
  transform = identity
)

write.csv(final_result, file = file.path(layer_output_dir, "platelet_age_gender_FDR.csv"))
