# Metabolite ~ Age association analysis (adjusted for Gender)
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

layer_output_dir <- file.path(project_dir, "output", "03_metabolite")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

metabolite <- read_tsv(file.path(raw_data_dir, "metabolic", "500FG raw metabolite_fromJianbo.tsv"))
metabolite_num <- as.data.frame(metabolite[, 22:479])
rownames(metabolite_num) <- metabolite$metabolite_identification
metabolite_trans <- as.data.frame(t(metabolite_num))  # transpose to samples x metabolites
colnames(metabolite_trans) <- rownames(metabolite_num)
rownames(metabolite_trans) <- colnames(metabolite_num)

basicPhenos <- read.csv(basic_phenos_path, row.names = 1)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(metabolite_trans, basicPhenos)
metabolite_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 2: age + metabolite using FDR, age not scaled + gender
final_result <- run_age_gender_lm(
  metabolite_filter,
  basicPhenos_filter_match,
  transform = function(x) log(x + 1)
)

write.csv(final_result, file = file.path(layer_output_dir, "metabolite_age_linear_gender_FDR.csv"))

write_sig_results(
  final_result,
  file.path(layer_output_dir, "metabolite_age_linear_gender_FDR_sig.csv")
)
