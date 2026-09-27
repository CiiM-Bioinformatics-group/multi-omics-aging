# Microbiome pathway abundance ~ Age association analysis (adjusted for
# Gender). Run from the repository root. No data is included in this repo
# -- see README for the expected input layout.
#
# NOTE (see README "Known issues"), kept as-is to match the original
# results:
# 1) if any pathway triggers the low-sample/zero-variance skip branch
#    inside run_age_gender_lm(), do.call(rbind, ...) errors on mismatched
#    columns -- a pre-existing issue in the original script.
# 2) the sig column is reclassified for the volcano plot colors before the
#    final sig-only export, so that export is always empty.
# 3) the volcano plot is built but never written to a file (no
#    png()/dev.off() around it, unlike the olink/hormone scripts).

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
library(car)

source("R/utils.R")

layer_output_dir <- file.path(project_dir, "output", "10_microbiome")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

microbiome <- read.table(
  file.path(raw_data_dir, "Molecular_phenotype", "500FG_microbiome_pathways.txt"),
  header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE
)

basicPhenos <- read.csv(basic_phenos_path)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(microbiome, basicPhenos)
microbiome_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 2: age + microbiome pathway using FDR, age not scaled + gender.
# asinh transform (raw values can be negative); skip features with fewer
# than 10 complete-case samples or zero variance after transform.
final_result <- run_age_gender_lm(
  microbiome_filter,
  basicPhenos_filter_match,
  transform = asinh,
  min_n = 10,
  require_variance = TRUE
)

write.csv(final_result, file = file.path(layer_output_dir, "microbiome_age_gender_FDR.csv"))
