# 3.1 Age prediction from each single omics layer (500FG, all samples of each layer)
#     Elastic net, 100 random 70/30 splits, 10-fold x 5 repeated CV for tuning.
#     DNA methylation is handled separately in 3.2-3.3.

source("prediction_utils.R")

# ---- Paths ----
data_dir   <- "data"
result_dir <- "results/prediction"
fig_dir    <- "results/figures"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pheno_age <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"))

# Replace Inf with the column's max finite value (used for cell counts)
replace_inf <- function(df) {
  out <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) col[is.infinite(col)] <- max(col[!is.infinite(col)], na.rm = TRUE)
    col
  }))
  rownames(out) <- rownames(df)
  out
}

# ---- Omics settings ----
omics_config <- list(
  cytokine = list(
    label = "Cytokine",
    read  = function() read.csv(file.path(data_dir, "cytokine_rank_normalized.csv"), row.names = 1),
    rds   = "cytokine_prediction_elastic_results.rds",
    pred  = "cytokine_pre_age_elasticnet_orign_sample_feature.csv"
  ),
  olink = list(
    label = "Proteomics",
    read  = function() read.csv(file.path(data_dir, "olink_rank_normalized.csv"), row.names = 1),
    pheno = function() {
      p <- read.csv(file.path(data_dir, "500FG_basicPhenos.csv"))
      p[!is.na(p$Age), ]
    },
    rds   = "olink_prediction_elastic_results.rds",
    pred  = "olink_pre_age_elasticnet_orign_sample_feature.csv"
  ),
  metabolite = list(
    label = "Metabolite",
    # samples in rows, metabolites in columns
    read  = function() read.csv(file.path(data_dir, "metabolite_rank_normalized.csv"),
                                row.names = 1, check.names = FALSE),
    rds   = "metabolite_prediction_results.rds",
    pred  = "metabolite_pre_age_elasticnet_orign_sample_feature.csv"
  ),
  cellcounts = list(
    label = "Cellcounts",
    read  = function() read.table(file.path(data_dir, "500FG_inverse_rank_normalized_cellcounts.txt"),
                                  header = TRUE, stringsAsFactors = FALSE),
    preprocess = replace_inf,
    rds   = "cellcount_prediction_elastic_results.rds",
    pred  = "cellcount_pre_age_elasticnet_orign_sample_feature.csv"
  ),
  microbiome = list(
    label = "Microbiome",
    read  = function() read.table(file.path(data_dir, "500FG_microbiome_pathways.txt"),
                                  header = TRUE, sep = "", stringsAsFactors = FALSE),
    rds   = "microbiome_prediction_elastic_results.rds",
    pred  = "microbiome_pre_age_elasticnet_orign_sample_feature.csv"
  )
)

# Choose which layers to run, e.g. omics_to_run <- "cytokine"
omics_to_run <- names(omics_config)

# ---- Run ----
for (name in omics_to_run) {
  cfg <- omics_config[[name]]
  message("Running: ", name)

  pheno <- if (is.null(cfg$pheno)) pheno_age else cfg$pheno()
  data  <- match_age(cfg$read(), pheno)
  if (!is.null(cfg$preprocess)) data <- cfg$preprocess(data)
  print(dim(data))

  results <- lapply(1:100, function(i) train_and_evaluate(data, seed = i, iteration = i))

  print(summarize_results(results))
  save_results(results, file.path(result_dir, cfg$rds), file.path(result_dir, cfg$pred))
  plot_r2(results, cfg$label, file.path(fig_dir, paste0("3.1_", name, "_cross_validation.png")))
}
