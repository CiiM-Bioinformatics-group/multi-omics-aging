# Cell counts (less-raw / inverse-rank-normalized subset) ~ Age association
# analysis (adjusted for Gender). 

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

layer_output_dir <- file.path(project_dir, "output", "07_cellcounts")
dir.create(layer_output_dir, showWarnings = FALSE, recursive = TRUE)

cellCounts <- read.csv(
  file.path(raw_data_dir, "Molecular_phenotype", "500FG_cellCounts_cellPerc.csv"),
  row.names = 1
)
cellCounts_select <- read.table(
  file.path(raw_data_dir, "Molecular_phenotype", "500FG_inverse_rank_normalized_cellcounts.txt"),
  header = TRUE, stringsAsFactors = FALSE
) # less-raw subset of cell types
cellCounts <- cellCounts[, colnames(cellCounts) %in% colnames(cellCounts_select)]

basicPhenos <- read.csv(basic_phenos_path)
rownames(basicPhenos) <- basicPhenos$ID_500fg

res <- prep_filter_match(cellCounts, basicPhenos)
cellCounts_filter <- res$feature_filter
basicPhenos_filter_match <- clean_gender(res$basicPhenos_filter_match)

# model 2: age + cell counts using FDR, age not scaled + gender
final_result <- run_age_gender_lm(
  cellCounts_filter,
  basicPhenos_filter_match,
  transform = function(x) log(x + 1)
)

# replace feature codes with their human-readable names
cellcounts_info <- read.table(file.path(raw_data_dir, "Molecular_phenotype", "500FG_cellcounts_info.txt"))
final_result$feature <- cellcounts_info$finalName[match(final_result$feature, cellcounts_info$code)]

write.csv(final_result, file = file.path(layer_output_dir, "cellcounts_age_rawless_gender_FDR.csv"))

write_sig_results(
  final_result,
  file.path(layer_output_dir, "cellcounts_age_rawless_gender_FDR_sig.csv")
)
