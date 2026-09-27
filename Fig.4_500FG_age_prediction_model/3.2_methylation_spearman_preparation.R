# 3.2 Spearman correlation of each CpG with age in the training set of each split
#     (all methylation samples, 100 splits). CpGs are processed in two halves to limit memory.
#     Output per half: list of 100 iterations with split indices and CpGs with p < 0.05.

library(caret)
library(dplyr)
library(parallel)

# ---- Paths ----
data_dir <- "data"
prep_dir <- "results/methylation_spearman"
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- 8

# ---- Load data ----
Mvalue_trans <- readRDS(file.path(data_dir, "Mvalue_trans.rds"))   # samples x CpGs
basicPhenos  <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"))

idx <- intersect(rownames(Mvalue_trans), basicPhenos$ID_500fg)
Mvalue_filter <- Mvalue_trans[rownames(Mvalue_trans) %in% idx, ]
age <- basicPhenos$Age[match(rownames(Mvalue_filter), basicPhenos$ID_500fg)]
rm(Mvalue_trans); gc()

# CpG halves (column ranges)
halves <- list(`1` = 1:427184, `2` = 427185:854368)

# ---- One split: Spearman test of every CpG vs age in the training set ----
spearman_split <- function(data, seed) {
  set.seed(seed)
  split_indices <- createDataPartition(data$Age, p = 0.7, list = FALSE)
  training <- data[split_indices, ]

  feature_names   <- setdiff(names(training), "Age")
  feature_batches <- split(feature_names, ceiling(seq_along(feature_names) / 500))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)

  spearman_results <- mclapply(feature_batches, function(batch) {
    sapply(batch, function(feature) {
      if (var(training[[feature]], na.rm = TRUE) == 0) return(c(cor = NA, pvalue = NA))
      test <- tryCatch(cor.test(training[[feature]], training$Age, method = "spearman"),
                       error = function(e) NULL)
      if (is.null(test)) return(c(cor = NA, pvalue = NA))
      c(cor = test$estimate, pvalue = test$p.value)
    }, simplify = FALSE)
  }, mc.cores = n_cores)

  spearman_pvalues <- unlist(lapply(spearman_results, function(batch) {
    sapply(batch, function(x) if (!is.null(x)) x["pvalue"] else NA)
  }))
  names(spearman_pvalues) <- unlist(lapply(spearman_results, names))

  selected_features <- names(spearman_pvalues[spearman_pvalues < 0.05])
  list(split_indices = split_indices,
       selected_features = selected_features,
       spearman_pvalues = spearman_pvalues[selected_features])
}

# ---- Run each half ----
for (h in names(halves)) {
  # Each half was originally run in a fresh R session: reset the RNG so splits are identical
  RNGkind("default", "default", "default")
  Mvalue_age <- cbind(Mvalue_filter[, halves[[h]]], Age = age)

  results <- lapply(1:100, function(i) spearman_split(Mvalue_age, seed = 123 + i))
  saveRDS(results, file.path(prep_dir, paste0("all_iterations_results", h, ".rds")))
  rm(Mvalue_age, results); gc()
}
