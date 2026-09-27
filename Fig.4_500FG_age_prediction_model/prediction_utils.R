# Shared functions for age prediction

library(caret)
library(dplyr)
library(glmnet)
library(ggplot2)

# Align an omics matrix (samples in rows) with phenotypes and append Age
match_age <- function(omics, pheno) {
  idx <- intersect(rownames(omics), pheno$ID_500fg)
  omics_filter <- omics[rownames(omics) %in% idx, ]
  pheno_match <- pheno[match(rownames(omics_filter), pheno$ID_500fg), ]
  stopifnot(identical(as.character(pheno_match$ID_500fg), rownames(omics_filter)))
  cbind(omics_filter, Age = pheno_match$Age)
}

# Fit elastic net on one train/validation split and evaluate
fit_and_evaluate <- function(training, validation, iteration,
                             alpha = seq(0.1, 0.9, by = 0.1),
                             lambda = seq(0, 1, by = 0.01),
                             keep_model = TRUE) {
  netGrid <- expand.grid(.alpha = alpha, .lambda = lambda)
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

  netFit <- train(Age ~ ., data = training, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")

  train_predictions <- predict(netFit, newdata = training)
  val_predictions   <- predict(netFit, newdata = validation)

  # Predicted age per sample
  pred_df <- rbind(
    data.frame(iteration = iteration, set = "train", sample_id = rownames(training),
               Age_true = training$Age, Age_pred = as.numeric(train_predictions)),
    data.frame(iteration = iteration, set = "validation", sample_id = rownames(validation),
               Age_true = validation$Age, Age_pred = as.numeric(val_predictions))
  )

  list(
    train_RMSE = sqrt(mean((training$Age - train_predictions)^2)),
    train_MSE  = mean((training$Age - train_predictions)^2),
    train_R2   = cor(training$Age, train_predictions)^2,
    val_RMSE   = sqrt(mean((validation$Age - val_predictions)^2)),
    val_MSE    = mean((validation$Age - val_predictions)^2),
    val_R2     = cor(validation$Age, val_predictions)^2,
    best_model = if (keep_model) netFit else NULL,
    predictions = pred_df
  )
}

# Random 70/30 split with a given seed, then fit
train_and_evaluate <- function(data, seed, iteration, alpha = seq(0.1, 0.9, by = 0.1)) {
  set.seed(seed)
  splitSample <- createDataPartition(data$Age, p = 0.7, list = FALSE)
  fit_and_evaluate(data[splitSample, ], data[-splitSample, ], iteration, alpha)
}

# Mean and SD of performance metrics across iterations
summarize_results <- function(results) {
  metrics <- c("train_RMSE", "val_RMSE", "train_MSE", "val_MSE", "train_R2", "val_R2")
  vals <- sapply(metrics, function(m) sapply(results, function(x) x[[m]]))
  data.frame(Metric = metrics, Mean = colMeans(vals), SD = apply(vals, 2, sd), row.names = NULL)
}

# Save results (.rds) and per-sample predictions (.csv)
save_results <- function(results, rds_file, pred_file) {
  saveRDS(results, rds_file)
  all_predictions <- do.call(rbind, lapply(results, function(x) x$predictions))
  write.csv(all_predictions, pred_file, row.names = FALSE)
}

# Boxplot of training vs validation R²
plot_r2 <- function(results, label, file) {
  plot_data_r2 <- data.frame(
    Set = rep(c("Training", "Validation"), each = length(results)),
    R2  = c(sapply(results, function(x) x$train_R2), sapply(results, function(x) x$val_R2))
  )
  g <- ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +
    labs(title = paste("R² Distribution Across Training and Validation in", label, "Sets"),
         x = label, y = "R²") +
    theme_minimal() +
    scale_fill_manual(values = c("lightblue", "lightgreen"))
  ggsave(file, g, width = 5, height = 4, dpi = 300, bg = "white")
  invisible(g)
}

# ---- Multi-omics helpers (3.4, 3.5) ----

# Six layers used in the multi-omics model (samples in rows)
load_multiomics_layers <- function(data_dir) {
  list(
    cytokine   = read.csv(file.path(data_dir, "cytokine_rank_normalized.csv"), row.names = 1),
    olink      = read.csv(file.path(data_dir, "olink_rank_normalized.csv"), row.names = 1),
    metabolite = read.csv(file.path(data_dir, "metabolite_rank_normalized.csv"), row.names = 1),
    cellcounts = read.csv(file.path(data_dir, "cellcounts_name_replace.csv"), row.names = 1),
    microbiome = read.csv(file.path(data_dir, "microbiome_rank_normalized.csv"), row.names = 1),
    Mvalue     = readRDS(file.path(data_dir, "Mvalue_trans.rds"))
  )
}

# Samples present in all layers and with phenotype, sorted
get_common_samples <- function(layers, pheno) {
  sort(Reduce(intersect, c(lapply(layers, rownames), list(pheno$ID_500fg))))
}

# Spearman correlation of every feature with Age in the training set, FDR (BH) across features
spearman_fdr_split <- function(data, split_indices, n_cores = 8) {
  training <- data[split_indices, , drop = FALSE]

  feature_names   <- setdiff(names(training), "Age")
  feature_batches <- split(feature_names, ceiling(seq_along(feature_names) / 500))

  spearman_results <- parallel::mclapply(feature_batches, function(batch) {
    sapply(batch, function(feature) {
      if (var(training[[feature]], na.rm = TRUE) == 0) return(c(cor = NA, pvalue = NA))
      test <- tryCatch(suppressWarnings(cor.test(training[[feature]], training$Age, method = "spearman")),
                       error = function(e) NULL)
      if (is.null(test)) return(c(cor = NA, pvalue = NA))
      c(cor = as.numeric(test$estimate), pvalue = test$p.value)
    }, simplify = FALSE)
  }, mc.cores = n_cores)

  feature_names <- unlist(lapply(spearman_results, names))
  spearman_pvalues <- unlist(lapply(spearman_results, function(b) sapply(b, function(x) x["pvalue"])))
  spearman_cor     <- unlist(lapply(spearman_results, function(b) sapply(b, function(x) x["cor"])))
  names(spearman_pvalues) <- names(spearman_cor) <- feature_names
  spearman_fdr <- p.adjust(spearman_pvalues, method = "fdr")

  list(split_indices = split_indices, feature_names = feature_names, spearman_cor = spearman_cor,
       spearman_pvalues = spearman_pvalues, spearman_fdr = spearman_fdr)
}
