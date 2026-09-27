# Immunoglobulin ~ Age association analysis (adjusted for Gender)


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

layer_output_dir <- file.path(project_dir, "output", "08_immunoglobulin")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

immunoglobulin <- read.table(
  file.path(raw_data_dir, "Molecular_phenotype", "500FG_log2_immunoglobulin_levels.txt"),
  header = TRUE, sep = "", stringsAsFactors = FALSE
)

basicPhenos <- read.csv(basic_phenos_path)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(immunoglobulin, basicPhenos)
immunoglobulin_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 2: age + immunoglobulin using FDR, age not scaled + gender
final_result <- run_age_gender_lm(
  immunoglobulin_filter,
  basicPhenos_filter_match,
  transform = identity
)

write.csv(final_result, file = file.path(layer_output_dir, "immunoglobulin_age_gender_FDR.csv"))

write_sig_results(
  final_result,
  file.path(layer_output_dir, "immunoglobulin_age_gender_FDR_sig.csv")
)
