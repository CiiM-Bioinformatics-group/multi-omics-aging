# Methylation (M-value) ~ Age association analysis (adjusted for Gender)


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

layer_output_dir <- file.path(project_dir, "output", "05_methylation")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

Mvalue <- readRDS(file.path(methylation_data_dir, "500FG.Mvalue.rds"))
Mvalue_trans <- as.data.frame(t(Mvalue))

basicPhenos <- read.csv(basic_phenos_path)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(Mvalue_trans, basicPhenos)
Mvalue_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 1: age + methylation (age not scaled) + Gender
final_result <- run_age_gender_lm(
  Mvalue_filter,
  basicPhenos_filter_match,
  transform = identity
)

saveRDS(final_result, file = file.path(layer_output_dir, "methylation_age_liner_gender_FDR.RDS"))
