# 3.4 Methylation feature selection for the multi-omics model (common samples of all six layers)
#     1) 100 fixed 70/30 splits of the common samples
#     2) Spearman correlation of each CpG with age in the training set, FDR within each CpG half
#     3) Keep CpGs with FDR < 0.05 (both halves merged) -> one JSON per iteration, used in 3.5
#     Common samples and splits are saved for 3.7.

source("prediction_utils.R")
library(parallel)
library(jsonlite)

# ---- Paths ----
data_dir <- "data"
fdr_dir  <- "results/fdr_common"
sel_dir  <- file.path(fdr_dir, "FDR_common005")
dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- 8
fdr_threshold <- 0.05

# ---- Common samples ----
layers <- load_multiomics_layers(data_dir)
basicPhenos <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"), row.names = 1)
common_samples <- get_common_samples(layers, basicPhenos)
age <- basicPhenos$Age[match(common_samples, basicPhenos$ID_500fg)]
Mvalue <- layers$Mvalue[common_samples, , drop = FALSE]
rm(layers); gc()
length(common_samples)
saveRDS(common_samples, file.path(fdr_dir, "common_samples.rds"))

# ---- 1) Fixed splits (shared by all layers of the multi-omics model) ----
RNGkind("L'Ecuyer-CMRG")
split_indices_list <- lapply(1:100, function(i) {
  set.seed(123 + i)
  as.integer(createDataPartition(age, p = 0.7, list = FALSE))
})
saveRDS(split_indices_list, file.path(fdr_dir, "split_indices_list.rds"))

# ---- 2) Spearman + FDR per split ----
# CpG halves (column ranges as in the original analysis; FDR is computed within each half)
halves <- list(`1` = 1:427184, `2` = 427184:854368)

for (h in names(halves)) {
  data_with_age <- as.data.frame(cbind(Mvalue[, halves[[h]]], Age = age))
  iteration_results <- lapply(1:100, function(i) spearman_fdr_split(data_with_age, split_indices_list[[i]], n_cores))
  write_json(iteration_results, file.path(fdr_dir, paste0("Methylation_FDR_common_results", h, ".json")), pretty = TRUE)
  rm(data_with_age, iteration_results); gc()
}

# ---- 3) Keep CpGs with FDR < 0.05 from both halves ----
results1 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results1.json"))
results2 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results2.json"))

feature_counts <- numeric(100)
for (i in 1:nrow(results1)) {
  fdr1 <- unlist(results1$spearman_fdr[[i]]); names1 <- unlist(results1$feature_names[[i]])
  fdr2 <- unlist(results2$spearman_fdr[[i]]); names2 <- unlist(results2$feature_names[[i]])
  sel1 <- names1[fdr1 < fdr_threshold]; sel_fdr1 <- fdr1[fdr1 < fdr_threshold]
  sel2 <- names2[fdr2 < fdr_threshold]; sel_fdr2 <- fdr2[fdr2 < fdr_threshold]

  merged_features <- unique(c(sel1, sel2))
  merged_fdr <- c(sel_fdr1, sel_fdr2)
  names(merged_fdr) <- c(sel1, sel2)
  merged_fdr <- tapply(merged_fdr, names(merged_fdr), min)[merged_features]   # smallest FDR per CpG
  feature_counts[i] <- length(merged_features)

  write_json(list(split_indices1 = results1$split_indices[[i]],
                  split_indices2 = results2$split_indices[[i]],
                  selected_features = merged_features,
                  fdr = as.numeric(merged_fdr)),
             file.path(sel_dir, sprintf("iteration_%d_fdr_info.json", i)), pretty = TRUE)
}
cat(sprintf("CpGs per iteration (FDR < %.2f): mean %.1f, range %d-%d\n",
            fdr_threshold, mean(feature_counts), min(feature_counts), max(feature_counts)))
