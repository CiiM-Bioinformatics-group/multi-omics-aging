## Age-explained variance by omics layer + sex (500FG aging)
##
## Stage 1 builds each layer's *_variance_input.csv from raw data
## Stage 2 uses those CSVs to produce:
##   - individual_layer_variance_gender_single.pdf
##   - cumulative_variance_gender_nature_order.pdf

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

## ---- paths (override via env vars; defaults are placeholders, not the
##      original cluster paths) ------------------------------------------
cfg <- list(
  input_dir       = Sys.getenv("VARIANCE_INPUT_DIR", "data/variance"),
  pheno_file      = Sys.getenv("VARIANCE_PHENO_FILE", "data/Age_group_basicPhenos.csv"),
  cellcounts_file = Sys.getenv("VARIANCE_CELLCOUNTS_FILE",
                                "data/500FG_inverse_rank_normalized_cellcounts.txt"),
  output_dir      = Sys.getenv("VARIANCE_OUTPUT_DIR", "output")
)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- function(name) file.path(cfg$input_dir, "cumulative_variance_input", name)

## ---- raw layer files (override via env vars; defaults are placeholders) ---
raw <- list(
  cytokine       = Sys.getenv("RAW_CYTOKINE_FILE", "data/raw/cytokine_levels.csv"),
  olink          = Sys.getenv("RAW_OLINK_FILE", "data/raw/olink_NPX.RDS"),
  metabolite     = Sys.getenv("RAW_METABOLITE_FILE", "data/raw/metabolite_levels.tsv"),
  hormone        = Sys.getenv("RAW_HORMONE_FILE", "data/raw/hormone_log2_levels.txt"),
  platelet       = Sys.getenv("RAW_PLATELET_FILE", "data/raw/platelet_log2_count.txt"),
  immunoglobulin = Sys.getenv("RAW_IMMUNOGLOBULIN_FILE", "data/raw/immunoglobulin_log2_levels.txt"),
  microbiome     = Sys.getenv("RAW_MICROBIOME_FILE", "data/raw/microbiome_pathways.txt")
)

## ---- Stage 1: raw data -> *_variance_input.csv ------------------------------
## Set VARIANCE_RUN_PREPROCESSING=FALSE to skip and reuse existing CSVs.
run_preprocessing <- toupper(Sys.getenv("VARIANCE_RUN_PREPROCESSING", "TRUE")) != "FALSE"

mean_impute <- function(x) apply(x, 2, function(col) ifelse(is.na(col), mean(col, na.rm = TRUE), col))
cap_inf     <- function(x) { x[is.infinite(x)] <- max(x[!is.infinite(x)], na.rm = TRUE); x }

preprocess_cytokine <- function() {
  d <- read.csv(raw$cytokine)
  rownames(d) <- d[, 1]
  d <- t(d[, -1])
  write.csv(mean_impute(d), input_file("cytokine_variance_input.csv"), row.names = TRUE)
}

preprocess_olink <- function() {
  d <- as.data.frame(readRDS(raw$olink))
  write.csv(mean_impute(d), input_file("olink_variance_input.csv"), row.names = TRUE)
}

preprocess_metabolite <- function() {
  d <- read_tsv(raw$metabolite, show_col_types = FALSE)
  d_num <- as.data.frame(d[, 22:479])
  rownames(d_num) <- d$metabolite_identification
  d_num <- t(d_num)
  write.csv(mean_impute(d_num), input_file("metabolite_variance_input.csv"), row.names = TRUE)
}

preprocess_hormone <- function() {
  d <- read.table(raw$hormone, header = TRUE, sep = "", stringsAsFactors = FALSE)
  write.csv(mean_impute(d), input_file("hormone_variance_input.csv"), row.names = TRUE)
}

preprocess_immunoglobulin <- function() {
  d <- read.table(raw$immunoglobulin, header = TRUE, sep = "", stringsAsFactors = FALSE)
  write.csv(cap_inf(mean_impute(d)), input_file("immunoglobulin_variance_input.csv"), row.names = TRUE)
}

## Drops columns that are >50% one repeated value or >50% missing, and
## zero-variance columns, before imputing (matches the ">50% missing
## value removed" note on this file's name).
preprocess_platelet <- function() {
  d <- read.table(raw$platelet, header = TRUE, sep = "", stringsAsFactors = FALSE)

  same_value <- sapply(d, function(x) max(table(x)) / length(x) > 0.5)
  d <- d[, !same_value, drop = FALSE]

  na_frac <- colSums(is.na(d)) / nrow(d)
  d <- d[, na_frac <= 0.5, drop = FALSE]

  low_var <- apply(d, 2, function(col) var(col, na.rm = TRUE) < 1e-6)
  d <- d[, !low_var, drop = FALSE]

  write.csv(cap_inf(mean_impute(d)), input_file("platelet_variance_input_removed.csv"), row.names = TRUE)
}

## Note: writes the raw table (not the imputed version), matching the
## source notebook.
preprocess_microbiome <- function() {
  d <- read.table(raw$microbiome, header = TRUE, sep = "", stringsAsFactors = FALSE)
  write.csv(d, input_file("microbiome_variance_input.csv"), row.names = TRUE)
}

if (run_preprocessing) {
  preprocess_cytokine()
  preprocess_olink()
  preprocess_metabolite()
  preprocess_hormone()
  preprocess_immunoglobulin()
  preprocess_platelet()
  preprocess_microbiome()
}

## ---- Stage 2: figures --------------------------------------------------

## ---- helpers ----------------------------------------------------------------
load_and_process <- function(file_path) {
  data <- read.csv(file_path)
  rownames(data) <- data[[1]]
  data[, -1, drop = FALSE]
}

clean_predictors <- function(data) {
  data <- as.data.frame(data)
  names(data) <- make.names(names(data), unique = TRUE)
  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  if (ncol(numeric_data) == 0) return(data.frame())
  numeric_data <- numeric_data[, colSums(!is.na(numeric_data)) > 0, drop = FALSE]
  if (ncol(numeric_data) == 0) return(data.frame())
  numeric_data[, apply(numeric_data, 2, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]
}

select_significant_features <- function(response, predictors, p_threshold = 0.05) {
  predictors <- clean_predictors(predictors)
  keep <- list()
  for (col in colnames(predictors)) {
    test <- tryCatch(cor.test(response, predictors[[col]], method = "spearman"),
                      error = function(e) NULL)
    if (!is.null(test) && !is.na(test$p.value) && test$p.value < p_threshold) {
      keep[[col]] <- predictors[[col]]
    }
  }
  if (length(keep) > 0) as.data.frame(keep) else data.frame()
}

remove_collinear_features <- function(predictors, threshold = 0.4) {
  predictors <- as.data.frame(predictors) %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(where(~ var(., na.rm = TRUE) > 0))
  if (ncol(predictors) < 2) return(predictors)

  cor_matrix <- cor(predictors, method = "spearman", use = "pairwise.complete.obs")
  to_remove <- c()
  for (i in 1:(ncol(cor_matrix) - 1)) {
    for (j in (i + 1):ncol(cor_matrix)) {
      if (abs(cor_matrix[i, j]) > threshold) {
        var_i <- var(predictors[[i]], na.rm = TRUE)
        var_j <- var(predictors[[j]], na.rm = TRUE)
        to_remove <- c(to_remove,
                        if (var_i < var_j) colnames(predictors)[i] else colnames(predictors)[j])
      }
    }
  }
  predictors[, !colnames(predictors) %in% unique(to_remove), drop = FALSE]
}

calculate_adjusted_r2 <- function(response, predictors = NULL, sex_df = NULL) {
  if (is.null(predictors)) predictors <- data.frame(row.names = seq_along(response))
  if (!is.data.frame(predictors)) predictors <- as.data.frame(predictors)

  model_df <- data.frame(response = response)
  if (!is.null(sex_df)) model_df <- cbind(model_df, sex_df)
  if (ncol(predictors) > 0) model_df <- cbind(model_df, predictors)

  model_df <- model_df[complete.cases(model_df), , drop = FALSE]
  colnames(model_df) <- make.names(colnames(model_df), unique = TRUE)
  if (ncol(model_df) <= 1) return(0)

  model <- lm(response ~ ., data = model_df)
  max(0, summary(model)$adj.r.squared)
}

relabel_layers <- function(x) {
  lookup <- c(
    "methylation" = "DNA methylation", "Proteomics" = "Proteomics",
    "Metabolites" = "Metabolomics", "Cytokine" = "Cytokine response",
    "Cell_count" = "Immune cell counts", "Microbiome" = "Microbiomes",
    "Hormones" = "Circulating endocrine traits",
    "Immunoglobulin" = "Circulating immune traits",
    "Platelets" = "Platelets", "Gender" = "Sex", "Sex" = "Sex"
  )
  out <- unname(lookup[x])
  ifelse(is.na(out), x, out)
}

## ---- data prep ----------------------------------------------------------------
load_all_layers <- function() {
  list(
    cytokine       = load_and_process(input_file("cytokine_variance_input.csv")),
    olink          = load_and_process(input_file("olink_variance_input.csv")),
    metabolite     = load_and_process(input_file("metabolite_variance_input.csv")),
    hormone        = load_and_process(input_file("hormone_variance_input.csv")),
    platelet       = load_and_process(input_file("platelet_variance_input_removed.csv")),
    immunoglobulin = load_and_process(input_file("immunoglobulin_variance_input.csv")),
    cellcounts     = read.table(cfg$cellcounts_file, header = TRUE, stringsAsFactors = FALSE),
    microbiome     = load_and_process(input_file("microbiome_variance_input.csv")),
    methylation    = { m <- load_and_process(input_file("Mvalue_pca_scores1.csv"))
                       colnames(m) <- paste0("m", colnames(m)); m }
  )
}

prepare_data <- function(layers, pheno_file) {
  basic_phenos <- load_and_process(pheno_file)

  all_data <- list(
    "methylation" = layers$methylation, "Proteomics" = layers$olink,
    "Cytokine" = layers$cytokine, "Metabolites" = layers$metabolite,
    "Cell_count" = layers$cellcounts, "Hormones" = layers$hormone,
    "Immunoglobulin" = layers$immunoglobulin, "Microbiome" = layers$microbiome,
    "Platelets" = layers$platelet
  )
  for (level in names(all_data)) {
    colnames(all_data[[level]]) <- make.unique(make.names(colnames(all_data[[level]])))
  }

  filter_hv <- function(data) data[grep("^HV", rownames(data)), , drop = FALSE]
  all_data <- lapply(all_data, filter_hv)

  common_samples <- Reduce(intersect, lapply(all_data, rownames))
  common_samples <- intersect(common_samples, rownames(basic_phenos))
  all_data <- lapply(all_data, function(x) x[common_samples, , drop = FALSE])
  basic <- basic_phenos[common_samples, , drop = FALSE]

  basic$Gender <- trimws(as.character(basic$Gender))
  basic$Gender[basic$Gender == ""] <- NA
  basic$Gender <- factor(basic$Gender)
  keep <- complete.cases(basic[, c("Age", "Gender")])
  basic <- basic[keep, , drop = FALSE]
  all_data <- lapply(all_data, function(x) x[rownames(basic), , drop = FALSE])

  list(all_data = all_data, basic = basic, response = basic$Age,
       sex_df = data.frame(Gender = basic$Gender))
}

## ---- scoring: cumulative (stacked) and individual (each layer + Sex) --------
compute_cumulative_r2 <- function(response, ordered_layers) {
  results <- data.frame()
  predictors_so_far <- data.frame(row.names = seq_along(response))
  cumulative_r2 <- calculate_adjusted_r2(response, predictors = NULL)
  results <- rbind(results, data.frame(Data_Level = "Baseline",
                                        Adjusted_R2 = cumulative_r2,
                                        Incremental_R2 = cumulative_r2))

  for (level in names(ordered_layers)) {
    layer_data <- ordered_layers[[level]]
    if (nrow(layer_data) == 0 || ncol(layer_data) == 0) next

    if (level == "Sex") {
      mm <- model.matrix(~ Sex, data = layer_data)
      predictors <- as.data.frame(mm[, -1, drop = FALSE], row.names = rownames(layer_data))
      if (ncol(predictors) == 0) next
    } else {
      predictors <- clean_predictors(layer_data)
      if (ncol(predictors) == 0) next
      predictors <- select_significant_features(response, predictors, p_threshold = 0.05)
      if (ncol(predictors) == 0) next
      predictors <- remove_collinear_features(predictors, threshold = 0.4)
      if (ncol(predictors) == 0) next
    }

    combined <- cbind(predictors_so_far, predictors)
    colnames(combined) <- make.unique(colnames(combined))

    adj_r2 <- calculate_adjusted_r2(response, combined)
    incremental_r2 <- max(0, adj_r2 - cumulative_r2)
    cumulative_r2 <- adj_r2

    results <- rbind(results, data.frame(Data_Level = level,
                                          Adjusted_R2 = cumulative_r2,
                                          Incremental_R2 = incremental_r2))
    predictors_so_far <- combined
  }
  results
}

compute_individual_r2 <- function(response, layers, sex_df) {
  out <- data.frame(Data_Level = "Sex",
                     Layer_Adjusted_R2 = calculate_adjusted_r2(response, predictors = NULL, sex_df = sex_df),
                     N_features = 0)

  for (level in names(layers)) {
    predictors <- clean_predictors(layers[[level]])
    if (ncol(predictors) == 0) next
    predictors <- select_significant_features(response, predictors, p_threshold = 0.05)
    if (ncol(predictors) == 0) next
    predictors <- remove_collinear_features(predictors, threshold = 0.4)
    if (ncol(predictors) == 0) next

    layer_r2 <- calculate_adjusted_r2(response, predictors, sex_df = sex_df)
    out <- rbind(out, data.frame(Data_Level = level,
                                  Layer_Adjusted_R2 = layer_r2,
                                  N_features = ncol(predictors)))
  }
  out
}

## ---- theme --------------------------------------------------------------------
theme_nature <- function(base_size = 5) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid      = element_blank(),
      plot.margin     = unit(c(1, 1, 1, 1), "mm"),
      axis.text       = element_text(size = base_size),
      axis.title      = element_text(size = base_size),
      legend.title    = element_text(size = base_size),
      legend.text     = element_text(size = base_size),
      legend.key.size = unit(2, "mm"),
      plot.title      = element_text(size = base_size)
    )
}

## ---- run ------------------------------------------------------------------
layers  <- load_all_layers()
prepped <- prepare_data(layers, cfg$pheno_file)

## Figure 1: individual_layer_variance_gender_single.pdf
individual_results <- compute_individual_r2(prepped$response, prepped$all_data, prepped$sex_df)

plot_df <- individual_results %>%
  mutate(Data_Level = relabel_layers(Data_Level), Adjusted_R2 = Layer_Adjusted_R2) %>%
  arrange(desc(Adjusted_R2)) %>%
  mutate(Data_Level = factor(Data_Level, levels = unique(Data_Level)))

color_palette_fig1 <- c(
  "DNA methylation" = "#8DD3C7", "Proteomics" = "#FB8072",
  "Metabolomics" = "#E9C46A", "Cytokine response" = "#80B1D3",
  "Immune cell counts" = "#BEBADA", "Microbiomes" = "#FCCDE5",
  "Circulating endocrine traits" = "#D98C6B",
  "Circulating immune traits" = "#5FAF9D",
  "Platelets" = "#D9D9D9", "Sex" = "#BC80BD"
)[levels(plot_df$Data_Level)]

p_individual <- ggplot(plot_df, aes(x = Data_Level, y = Adjusted_R2, fill = Data_Level)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_fill_manual(values = color_palette_fig1) +
  theme_nature(base_size = 5) +
  theme(
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none",
    aspect.ratio = 0.7
  ) +
  labs(y = expression("Individual adjusted R"^2))

ggsave(file.path(cfg$output_dir, "Fig2b_individual_layer_variance_gender_single.pdf"),
       plot = p_individual, width = 80, height = 61, units = "mm")

## Figure 2: cumulative_variance_gender_order.pdf
ordered_levels <- c("methylation", "Proteomics", "Metabolites", "Cytokine",
                     "Hormones", "Immunoglobulin", "Cell_count", "Sex",
                     "Microbiome", "Platelets")

cumulative_layers <- prepped$all_data
cumulative_layers$Sex <- data.frame(Sex = prepped$basic$Gender, row.names = rownames(prepped$basic))
cumulative_layers <- cumulative_layers[ordered_levels]

results <- compute_cumulative_r2(prepped$response, cumulative_layers) %>%
  dplyr::filter(Data_Level != "Baseline")

results <- results %>%
  mutate(Data_Level = factor(Data_Level, levels = unique(Data_Level)),
         Age = "Age Category",
         Data_Level = relabel_layers(as.character(Data_Level)),
         Data_Level = factor(Data_Level, levels = unique(Data_Level)))

color_palette_fig2 <- c(
  "Sex" = "#BC80BD", "DNA methylation" = "#8DD3C7", "Proteomics" = "#FB8072",
  "Metabolomics" = "#E9C46A", "Cytokine response" = "#80B1D3",
  "Immune cell counts" = "#BEBADA", "Microbiomes" = "#FCCDE5",
  "Circulating endocrine traits" = "#D98C6B",
  "Circulating immune traits" = "#5FAF9D"
)[levels(results$Data_Level)]

results$Data_Level <- factor(results$Data_Level, levels = rev(levels(results$Data_Level)))
color_palette_fig2 <- color_palette_fig2[levels(results$Data_Level)]

p_cumulative <- ggplot(results, aes(x = Age, y = Incremental_R2, fill = Data_Level)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE), width = 0.6) +
  scale_fill_manual(values = color_palette_fig2) +
  scale_x_discrete(labels = c("Age Category" = "Age")) +
  theme_nature(base_size = 5) +
  theme(
    axis.text.x = element_text(size = 5, face = "bold"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 3
  ) +
  labs(y = "Incremental Adjusted R²", fill = "Category added")

ggsave(file.path(cfg$output_dir, "Fig2c_cumulative_variance_gender_order.pdf"),
       plot = p_cumulative, width = 70, height = 88 * 7 / 10, units = "mm")
