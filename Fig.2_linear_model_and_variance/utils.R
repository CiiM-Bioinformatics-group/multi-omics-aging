# Path config + shared helpers for the 1.x_*.R layer scripts.
# Run scripts from the repository root: source("R/utils.R")

## --- paths (edit defaults, or override via env vars) ---
# raw_data_dir: shared 500FG cohort raw-data store (read-only)
# project_dir : this project's own derived data + outputs
raw_data_dir <- Sys.getenv("FG500_RAW_DATA_DIR", unset = "data/raw_cohort")
project_dir  <- Sys.getenv("FG500_PROJECT_DIR",  unset = ".")

# Age/Gender/ID phenotype table shared by every layer
basic_phenos_path <- file.path(project_dir,"Age_group_basicPhenos.csv")

## --- shared functions ---

# Restrict feature_mat + basicPhenos to shared sample IDs, aligned to
# feature_mat's row order.
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- intersect(rownames(feature_mat), rownames(basicPhenos))
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}

# Trim whitespace in Gender; blank -> NA.
clean_gender <- function(basicPhenos_filter_match) {
  basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
  basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA
  basicPhenos_filter_match
}

# Fit feature ~ age + Gender per column of feature_filter; return
# estimate/p/signed -log10(p)/padj (BH)/sig ("sig" if padj < 0.05).
#
# transform: applied to each feature column before modeling -- differs by
# layer (log(x+1), identity, or asinh) and changes the numeric results, so
# match it to the original per-layer script.
# min_n / require_variance: microbiome-only. Drop incomplete rows first;
# skip (NA row) a feature with fewer than min_n samples or zero variance.
run_age_gender_lm <- function(feature_filter,
                               basicPhenos_filter_match,
                               transform = function(x) log(x + 1),
                               min_n = 1,
                               require_variance = FALSE) {

  all_results <- vector("list", ncol(feature_filter))

  for (i in seq_len(ncol(feature_filter))) {
    data <- data.frame(
      feature = transform(feature_filter[, i]),
      age = basicPhenos_filter_match$Age,
      Gender = factor(basicPhenos_filter_match$Gender)
    )

    if (require_variance || min_n > 1) {
      data <- data[complete.cases(data), ]
      if (nrow(data) < min_n || sd(data$feature) == 0) {
        all_results[[i]] <- data.frame(
          feature  = colnames(feature_filter)[i],
          estimate = NA_real_,
          p        = NA_real_
        )
        next
      }
    }

    mod <- lm(feature ~ age + Gender, data = data)
    cf <- summary(mod)$coefficients
    res <- cf["age", c("Estimate", "Pr(>|t|)")]

    result <- data.frame(
      feature  = colnames(feature_filter)[i],
      estimate = res["Estimate"],
      p        = res["Pr(>|t|)"]
    )
    result$value <- -log10(result$p) * result$estimate

    all_results[[i]] <- result
  }

  final_result <- do.call(rbind, all_results)
  final_result$padj <- p.adjust(final_result$p, method = "BH")
  final_result$sig <- ifelse(final_result$padj < 0.05, "sig", "no")
  final_result
}

# Filter to sig rows (feature, estimate, p, padj) and write them out.
write_sig_results <- function(final_result, path) {
  final_result_sig <- final_result %>%
    dplyr::filter(sig == "sig") %>%
    dplyr::select(feature, estimate, p, padj)
  readr::write_excel_csv(final_result_sig, path)
  final_result_sig
}
