# 3.7 Spearman + FDR feature ranking for the non-methylation layers (common samples, splits from 3.4)
#     Output: {Layer}_FDR_common_results.json (used in 3.8), plus CSV copies of inputs for the Python scripts.

source("prediction_utils.R")
library(jsonlite)

# ---- Paths ----
data_dir <- "data"
fdr_dir  <- "results/fdr_common"
n_cores  <- 8

# ---- Load data ----
cytokine    <- readRDS(file.path(data_dir, "cytokine_filter_transposed.rds"))
olink       <- readRDS(file.path(data_dir, "NPX_FG500.RDS"))
metabolite  <- readRDS(file.path(data_dir, "metabolite_trans.rds"))
cellcounts  <- read.csv(file.path(data_dir, "cellcounts_name_replace.csv"), row.names = 1)
microbiome  <- read.table(file.path(data_dir, "500FG_microbiome_pathways.txt"),
                          header = TRUE, sep = "", stringsAsFactors = FALSE)
basicPhenos <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"), row.names = 1)
rownames(basicPhenos) <- basicPhenos$ID_500fg

# Common samples and fixed splits from 3.4
common_samples     <- readRDS(file.path(fdr_dir, "common_samples.rds"))
split_indices_list <- readRDS(file.path(fdr_dir, "split_indices_list.rds"))

filter_samples <- function(data, samples) {
  data <- data[rownames(data) %in% samples, , drop = FALSE]
  data[samples, , drop = FALSE]
}
basicPhenos_filter <- filter_samples(basicPhenos, common_samples)

datasets_to_process <- list(
  Cytokine   = filter_samples(cytokine, common_samples),
  Proteomics = filter_samples(olink, common_samples),
  Metabolite = filter_samples(metabolite, common_samples),
  Microbiome = filter_samples(microbiome, common_samples),
  Cellcounts = filter_samples(cellcounts, common_samples)
)

# ---- Spearman + FDR per layer and split ----
for (layer_name in names(datasets_to_process)) {
  data <- datasets_to_process[[layer_name]]
  stopifnot(all(rownames(data) == rownames(basicPhenos_filter)))
  data_with_age <- as.data.frame(cbind(data, Age = basicPhenos_filter$Age))

  iteration_results <- lapply(1:100, function(i) spearman_fdr_split(data_with_age, split_indices_list[[i]], n_cores))
  write_json(iteration_results, file.path(fdr_dir, paste0(layer_name, "_FDR_common_results.json")), pretty = TRUE)
}

# ---- CSV copies of the inputs read by the Python scripts (3.9, 3.12) ----
write.csv(cytokine,   file.path(data_dir, "cytokine_filter_transposed.csv"))
write.csv(olink,      file.path(data_dir, "olink_NPX_FG500.csv"))
write.csv(metabolite, file.path(data_dir, "metabolite_trans.csv"))
write.csv(readRDS(file.path(data_dir, "Mvalue_trans.rds")), file.path(data_dir, "Mvalue_trans.csv"))
