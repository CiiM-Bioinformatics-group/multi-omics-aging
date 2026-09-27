# 3.3 Age prediction from genome-wide DNA methylation (all methylation samples)
#     CpGs pre-selected per split by Spearman correlation with age (p < 1e-4, from 3.2), then elastic net.

source("prediction_utils.R")

# ---- Paths ----
data_dir   <- "data"
result_dir <- "results/prediction"
prep_dir   <- "results/methylation_spearman"
fig_dir    <- "results/figures"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
Mvalue      <- readRDS(file.path(data_dir, "Mvalue_trans.rds"))   # samples x CpGs
spearman1   <- readRDS(file.path(prep_dir, "all_iterations_results1.rds"))
spearman2   <- readRDS(file.path(prep_dir, "all_iterations_results2.rds"))
basicPhenos <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"))

Mvalue_age <- match_age(Mvalue, basicPhenos)
rm(Mvalue); gc()
dim(Mvalue_age)

# ---- One split: Spearman-selected CpGs + elastic net ----
train_and_evaluate_spearman <- function(data, seed, iteration, pvalue_threshold) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  split_indices <- createDataPartition(data$Age, p = 0.7, list = FALSE)

  # Spearman results for this iteration (CpGs were tested in two batches)
  s1 <- spearman1[[iteration]]
  s2 <- spearman2[[iteration]]
  all_pvalues  <- c(s1$spearman_pvalues, s2$spearman_pvalues)
  all_features <- c(s1$selected_features, s2$selected_features)

  # Use the same split as the Spearman pre-selection
  if (!identical(split_indices, s1$split_indices)) {
    warning("Split indices mismatch. Using split indices from Spearman results.")
    split_indices <- s1$split_indices
  }

  selected_features <- all_features[all_pvalues < pvalue_threshold]
  if (length(selected_features) == 0) {
    stop(paste("No features selected at p <", pvalue_threshold, "in iteration", iteration))
  }

  cols <- c(selected_features, "Age")
  res <- fit_and_evaluate(data[split_indices, cols, drop = FALSE],
                          data[-split_indices, cols, drop = FALSE],
                          iteration, keep_model = FALSE)
  gc()
  c(res, list(selected_features = selected_features))
}

# ---- Run 100 iterations ----
results <- lapply(1:100, function(i) {
  train_and_evaluate_spearman(Mvalue_age, seed = 123 + i, iteration = i, pvalue_threshold = 1e-04)
})

# ---- Summary ----
num_selected_features <- sapply(results, function(x) length(x$selected_features))
print(rbind(summarize_results(results),
            data.frame(Metric = "n_selected_features",
                       Mean = mean(num_selected_features), SD = sd(num_selected_features))))

# ---- Save and plot ----
save_results(results,
             file.path(result_dir, "methylation_prediction_pvalue_select1e04_results_spearman_elasticnet.rds"),
             file.path(result_dir, "methylation_spearman_pre_age_elasticnet.csv"))
plot_r2(results, "Methylation", file.path(fig_dir, "3.3_methylation_spearman_cross_validation.png"))
