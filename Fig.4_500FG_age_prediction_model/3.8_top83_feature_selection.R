# 3.8 Top 83 features per layer (smallest Spearman FDR) for each split, used by the TabPFN models in 3.9
#     Non-methylation layers from 3.7; methylation from the two CpG halves of 3.4.

library(jsonlite)

# ---- Paths ----
fdr_dir <- "results/fdr_common"
top_dir <- file.path(fdr_dir, "FDR_top83")
n_top   <- 83

# Keep the n_top features with the smallest FDR and write one JSON per iteration
write_top <- function(fdr, feature_names, split_indices, out_dir, i) {
  order_idx <- order(fdr)
  n_select  <- min(n_top, length(feature_names))
  result <- list(split_indices = split_indices,
                 selected_features = feature_names[order_idx][1:n_select],
                 fdr = fdr[order_idx][1:n_select])
  write_json(result, file.path(out_dir, sprintf("iteration_%d_fdr_info.json", i)), pretty = TRUE)
}

# ---- Non-methylation layers ----
for (layer_name in c("Cytokine", "Proteomics", "Metabolite", "Microbiome", "Cellcounts")) {
  all_results <- fromJSON(file.path(fdr_dir, paste0(layer_name, "_FDR_common_results.json")))
  out_dir <- file.path(top_dir, paste0(layer_name, "_spearman_preparation_top83"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in 1:nrow(all_results)) {
    write_top(unlist(all_results$spearman_fdr[[i]]), unlist(all_results$feature_names[[i]]),
              all_results$split_indices[[i]], out_dir, i)
  }
}

# ---- Methylation (both CpG halves pooled) ----
all_results1 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results1.json"))
all_results2 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results2.json"))
out_dir <- file.path(top_dir, "Methylation_merged_spearman_preparation_top83")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (i in 1:min(nrow(all_results1), nrow(all_results2))) {
  write_top(c(unlist(all_results1$spearman_fdr[[i]]), unlist(all_results2$spearman_fdr[[i]])),
            c(unlist(all_results1$feature_names[[i]]), unlist(all_results2$feature_names[[i]])),
            all_results1$split_indices[[i]], out_dir, i)
}
