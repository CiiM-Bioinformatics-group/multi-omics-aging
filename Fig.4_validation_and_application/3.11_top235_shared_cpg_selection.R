# 3.11 Methylation features for the 300BCG validation: among CpGs shared with 300BCG (3.10),
#      FDR recomputed from the 500FG Spearman p-values (3.4) and the top 235 CpGs kept per split.
#      Run from the repository root: Rscript validation/<this script>

library(jsonlite)

# ---- Paths ----
fdr_dir <- "results/fdr_common"
bcg_dir <- "data/300bcg"
out_dir <- "results/fdr_shared_top235/Methylation_FDR_shared_top235"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n_top <- 235

# ---- Load data ----
common_cols  <- read.csv(file.path(bcg_dir, "common_cols.csv"))
all_results1 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results1.json"))
all_results2 <- fromJSON(file.path(fdr_dir, "Methylation_FDR_common_results2.json"))

# ---- Top 235 shared CpGs per split ----
feature_counts <- numeric(100)
for (i in 1:min(nrow(all_results1), nrow(all_results2))) {
  spearman_pvalues <- c(unlist(all_results1$spearman_pvalues[[i]]), unlist(all_results2$spearman_pvalues[[i]]))
  feature_names    <- c(unlist(all_results1$feature_names[[i]]), unlist(all_results2$feature_names[[i]]))

  # Shared CpGs only, then FDR within this set
  keep <- substr(feature_names, 1, 10) %in% common_cols$x
  filtered_names <- feature_names[keep]
  filtered_fdr   <- p.adjust(as.numeric(spearman_pvalues[keep]), method = "fdr")

  order_idx <- order(filtered_fdr)
  n_select  <- min(n_top, length(filtered_names))
  result <- list(split_indices = all_results1$split_indices[[i]],
                 selected_features = filtered_names[order_idx][1:n_select],
                 fdr = filtered_fdr[order_idx][1:n_select])
  write_json(result, file.path(out_dir, sprintf("iteration_%d_fdr_info.json", i)), pretty = TRUE)
  feature_counts[i] <- n_select
}
cat("CpGs per iteration:", mean(feature_counts), "\n")
