# 3.5 Multi-omics age prediction (common samples of all six layers)
#     Cytokine, proteomics, metabolomics, cell counts and microbiome: all features;
#     methylation: CpGs with FDR < 0.05 per split (from 3.4, top 10,000 by FDR if > 20,000).
#     Elastic net on the same 100 fixed splits as 3.4.

source("prediction_utils.R")
library(jsonlite)

# ---- Paths ----
data_dir   <- "data"
sel_dir    <- "results/fdr_common/FDR_common005"
result_dir <- "results/prediction"
fig_dir    <- "results/figures"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data and build the combined table ----
layers <- load_multiomics_layers(data_dir)
basicPhenos <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"), row.names = 1)
common_samples <- get_common_samples(layers, basicPhenos)

combined_table <- do.call(cbind, unname(lapply(layers, function(x) x[common_samples, , drop = FALSE])))
rownames(combined_table) <- common_samples
combined_age <- cbind(combined_table, Age = basicPhenos$Age[match(common_samples, basicPhenos$ID_500fg)])
rm(combined_table); gc()
dim(combined_age)

non_methyl_features <- unlist(lapply(layers[names(layers) != "Mvalue"], colnames))
layer_features <- lapply(layers[names(layers) != "Mvalue"], colnames)
rm(layers); gc()

# ---- 100 iterations ----
results <- lapply(1:100, function(i) {
  Mvalue_info   <- fromJSON(file.path(sel_dir, sprintf("iteration_%d_fdr_info.json", i)))
  split_indices <- as.integer(unlist(Mvalue_info$split_indices1))

  # Methylation: FDR-selected CpGs; if > 20,000, keep the 10,000 with the smallest FDR
  Mvalue_features <- Mvalue_info$selected_features
  if (length(Mvalue_features) > 20000) {
    names(Mvalue_info$fdr) <- Mvalue_info$selected_features
    Mvalue_features <- names(sort(Mvalue_info$fdr, decreasing = FALSE))[1:10000]
  }

  combined_features  <- setdiff(unique(c(non_methyl_features, Mvalue_features)), "Age")
  available_features <- intersect(combined_features, colnames(combined_age))
  features_count <- c(sapply(layer_features, function(f) length(intersect(f, available_features))),
                      Mvalue = length(intersect(Mvalue_features, available_features)))

  cols <- c(available_features, "Age")
  set.seed(123 + i)   # fixes the CV folds used for tuning
  res <- fit_and_evaluate(combined_age[split_indices, cols, drop = FALSE],
                          combined_age[-split_indices, cols, drop = FALSE],
                          iteration = i, lambda = 10^seq(-4, 1, length.out = 100))

  coefficients <- as.matrix(coef(res$best_model$finalModel, s = res$best_model$bestTune$lambda))
  selected <- rownames(coefficients)[coefficients[, 1] != 0]
  c(res, list(coefficients = coefficients,
              selected_features = setdiff(selected, "(Intercept)"),
              features_count = features_count))
})

# ---- Summary ----
n_features <- sapply(results, function(x) sum(x$features_count))
print(summarize_results(results))
cat("Input features per iteration: mean", mean(n_features), "min", min(n_features), "max", max(n_features), "\n")

# ---- Save and plot ----
save_results(results,
             file.path(result_dir, "multi_prediction_results_original_methy_separ_elasticnet.rds"),
             file.path(result_dir, "multi_prediction_age_original_methy_separ_elasticnet.csv"))
plot_r2(results, "Multi-omics", file.path(fig_dir, "3.5_multiomics_cross_validation.png"))
