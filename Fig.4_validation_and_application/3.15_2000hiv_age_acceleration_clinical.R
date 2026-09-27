# 3.15 2000HIV: TabPFN age acceleration and its association with clinical variables
#      Age acceleration = residual of mean predicted age (100 TabPFN models, 3.14) on chronological age;
#      compared with epigenetic-clock age accelerations (Horvath, Hannum, PhenoAge, GrimAge).
#      Run from the repository root: Rscript validation/3.15_2000hiv_age_acceleration_clinical.R

library(dplyr)
library(ggplot2)

# ---- Paths ----
hiv_dir <- "data/2000hiv"
out_dir <- "results/validation_2000hiv"

# ---- TabPFN predicted age and age acceleration ----
pred_all <- read.csv(file.path(out_dir, "HIV_all_models_all_iterations_predicted_age.csv"))
tabpfn_mean <- pred_all[pred_all$Model == "TabPFN", ] %>%
  group_by(sample_id) %>%
  summarise(Predicted_Age = mean(predicted_age, na.rm = TRUE), Actual_Age = first(true_age), .groups = "drop")
tabpfn_mean <- tabpfn_mean[complete.cases(tabpfn_mean[, c("Predicted_Age", "Actual_Age")]), ]
tabpfn_mean$Age_Acceleration <- residuals(lm(Predicted_Age ~ Actual_Age, data = tabpfn_mean))

tabpfn_AA <- data.frame(Sample_ID = tabpfn_mean$sample_id, Predicted_Age = tabpfn_mean$Predicted_Age,
                        Actual_Age = tabpfn_mean$Actual_Age, Age_Acceleration = tabpfn_mean$Age_Acceleration)
rownames(tabpfn_AA) <- tabpfn_AA$Sample_ID

cat(sprintf("TabPFN in 2000HIV: R2=%.4f | RMSE=%.2f | MAE=%.2f | N=%d\n",
            cor(tabpfn_AA$Actual_Age, tabpfn_AA$Predicted_Age)^2,
            sqrt(mean((tabpfn_AA$Actual_Age - tabpfn_AA$Predicted_Age)^2)),
            mean(abs(tabpfn_AA$Actual_Age - tabpfn_AA$Predicted_Age)), nrow(tabpfn_AA)))
write.csv(tabpfn_AA, file.path(out_dir, "full_model_tabpfn_predicted_age_2000hiv_residual.csv"))

# ---- Phenotypes and epigenetic clocks ----
phenotype <- read.csv2(file.path(hiv_dir, "220317_2000hiv_study_export_processed.csv"))
rownames(phenotype) <- phenotype$Record.Id

AA <- read.csv(file.path(hiv_dir, "Aging.200HIV.csv"), row.names = 1)
AA_filter <- AA %>% select(EAA_Horvath = AgeAccelerationResidual, EAA_Hannum = AgeAccelerationResidualHannum,
                           EAA_Phenoage = AgeAccelPheno, EAA_Grimage = AgeAccelGrim)

cs <- Reduce(intersect, list(rownames(tabpfn_AA), rownames(AA_filter), rownames(phenotype)))
pheno_raw_match <- phenotype[cs, , drop = FALSE]
names(pheno_raw_match) <- trimws(names(pheno_raw_match))

# ---- Clinical variables ----
required_raw_columns <- c(
  "Record.Id", "AGE", "SEX_BIRTH", "Institute.Abbreviation", "DATE_VISIT", "SMOKING", "SMOKING_PY",
  "BMI_BASELINE", "RRSYST1", "RRSYST2", "RRSYST3", "RRDIA1", "RRDIA2", "RRDIA3", "LAB_eGFR",
  "FIB_CAP_MED", "MH_TR_END_D.DM1", "MH_TR_END_D.DM2", "MH_TR_CV_D.Hypertensie", "HYPERTENSION_YEAR",
  "MH_TR_CV_D.Myocard_infarct", "MH_TR_CV_D.Stroke", "MH_TR_CV_D.PAV", "MH_TR_CV_D.Angina_pectoris",
  "HIV_DATE", "HIV_DIAG", "DATE_START_CART", "YEAR_START_CART",
  "CD4_NADIR", "CD4_LATEST", "CD8_LATEST", "VL_UND_DATE", "VL_LATEST_DATE", "RESIDUAL_VIR", "CART_MONODUO"
)

# Alternative column names used in the export
alias_map <- list(
  "MH_TR_CV_D.Hypertensie"     = c("MH_TR_CV_D.Hypertension", "CVDRISC_Hypertensie"),
  "MH_TR_CV_D.Myocard_infarct" = c("MH_TR_CV_D.Myocardial_infarction", "MH_TR_CV_D.Myocard_infarction")
)
for (target in names(alias_map)) {
  if (!(target %in% names(pheno_raw_match))) {
    hit <- alias_map[[target]][alias_map[[target]] %in% names(pheno_raw_match)]
    if (length(hit) > 0) pheno_raw_match[[target]] <- pheno_raw_match[[hit[1]]]
  }
}
for (target in required_raw_columns) {
  if (!(target %in% names(pheno_raw_match))) {
    idx <- which(tolower(names(pheno_raw_match)) == tolower(target))
    if (length(idx) > 0) pheno_raw_match[[target]] <- pheno_raw_match[[idx[1]]]
  }
}
miss <- setdiff(required_raw_columns, names(pheno_raw_match))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
to_bin <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "yes", "y", "true", "positive", "present")] <- 1L
  out[z %in% c("0", "no", "n", "false", "negative", "absent")] <- 0L
  out
}
to_date <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA
  as.Date(x_chr, tryFormats = c("%Y-%m-%d", "%d-%m-%Y", "%Y/%m/%d", "%d/%m/%Y",
                                "%Y-%m-%d %H:%M:%S", "%d-%m-%Y %H:%M:%S"))
}
row_mean_na <- function(...) {
  m <- cbind(...)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  out <- rowMeans(m, na.rm = TRUE)
  out[rowSums(!is.na(m)) == 0] <- NA_real_
  out
}
row_any_yes <- function(...) {
  m <- cbind(...)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  obs <- rowSums(!is.na(m)) > 0
  out <- rep(NA_integer_, nrow(m))
  out[obs] <- as.integer(rowSums(m[obs, , drop = FALSE] == 1L, na.rm = TRUE) > 0)
  out
}
safe_year <- function(x) {
  y <- suppressWarnings(as.numeric(format(x, "%Y")))
  y[!is.finite(y)] <- NA_real_
  y
}

p <- pheno_raw_match
visit_date <- to_date(p$DATE_VISIT)
hiv_date   <- to_date(p$HIV_DATE)
art_date   <- to_date(p$DATE_START_CART)
vl_und_date <- to_date(p$VL_UND_DATE)

clinical_hiv_prepared <- tibble(
  participant_id = as.character(p$Record.Id),
  sex_at_birth = factor(ifelse(to_num(p$SEX_BIRTH) == 1, "Female", ifelse(to_num(p$SEX_BIRTH) == 0, "Male", NA)),
                        levels = c("Male", "Female")),
  smoking_status = factor(c("Never", "Current", "Former")[to_num(p$SMOKING) + 1], levels = c("Never", "Former", "Current")),
  BMI = {x <- to_num(p$BMI_BASELINE); x[x <= 0] <- NA; x},
  mean_sbp_mmHg = row_mean_na(to_num(p$RRSYST1), to_num(p$RRSYST2), to_num(p$RRSYST3)),
  mean_dbp_mmHg = row_mean_na(to_num(p$RRDIA1), to_num(p$RRDIA2), to_num(p$RRDIA3)),
  reduced_egfr_lt90 = ifelse(to_num(p$LAB_eGFR) == 0, 1L, ifelse(to_num(p$LAB_eGFR) == 1, 0L, NA)),
  fibroscan_cap_median = {x <- to_num(p$FIB_CAP_MED); x[x <= 0] <- NA; x},
  confirmed_diabetes = row_any_yes(to_bin(p$MH_TR_END_D.DM1), to_bin(p$MH_TR_END_D.DM2)),
  established_ascvd = row_any_yes(to_bin(p$MH_TR_CV_D.Myocard_infarct), to_bin(p$MH_TR_CV_D.Stroke),
                                  to_bin(p$MH_TR_CV_D.PAV), to_bin(p$MH_TR_CV_D.Angina_pectoris)),
  hiv_duration_years = {
    d1 <- as.numeric(visit_date - hiv_date) / 365.25
    d2 <- safe_year(visit_date) - to_num(p$HIV_DIAG)
    out <- ifelse(!is.na(d1), d1, d2); out[out < 0 | out > 100] <- NA; out
  },
  art_duration_years = {
    d1 <- as.numeric(visit_date - art_date) / 365.25
    d2 <- safe_year(visit_date) - to_num(p$YEAR_START_CART)
    out <- ifelse(!is.na(d1), d1, d2); out[out < 0 | out > 100] <- NA; out
  },
  cd4_nadir_cells_uL  = {x <- to_num(p$CD4_NADIR) * 1000; x[x < 0] <- NA; x},
  cd4_latest_cells_uL = {x <- to_num(p$CD4_LATEST) * 1000; x[x < 0] <- NA; x},
  cd8_latest = to_num(p$CD8_LATEST),
  CD4_CD8_ratio = {r <- to_num(p$CD4_LATEST) / to_num(p$CD8_LATEST); r[!is.finite(r)] <- NA_real_; r},
  consistent_unquantifiable_vl_3y = ifelse(to_num(p$RESIDUAL_VIR) == 1, 1L, ifelse(to_num(p$RESIDUAL_VIR) == 0, 0L, NA)),
  prior_mono_dual_art = to_bin(p$CART_MONODUO)
) %>%
  mutate(
    cd4_nadir_per100  = cd4_nadir_cells_uL / 100,
    cd4_latest_per100 = cd4_latest_cells_uL / 100,
    cd8_latest_per100 = cd8_latest * 10,   # raw CD8 in cells/nL
    across(c(reduced_egfr_lt90, confirmed_diabetes, established_ascvd,
             consistent_unquantifiable_vl_3y, prior_mono_dual_art), ~ factor(.x, levels = c(0, 1)))
  )

predictors <- c("sex_at_birth", "smoking_status", "BMI", "mean_sbp_mmHg", "mean_dbp_mmHg",
                "reduced_egfr_lt90", "fibroscan_cap_median", "confirmed_diabetes", "established_ascvd",
                "hiv_duration_years", "art_duration_years", "cd4_nadir_per100", "cd4_latest_per100",
                "cd8_latest_per100", "CD4_CD8_ratio", "consistent_unquantifiable_vl_3y", "prior_mono_dual_art")

# ---- Merge age accelerations with clinical variables ----
common_ids <- Reduce(intersect, list(rownames(AA_filter), rownames(tabpfn_AA), clinical_hiv_prepared$participant_id))
AA_match <- AA_filter[common_ids, , drop = FALSE]
AA_match$Multi_omics_AA <- to_num(tabpfn_AA[common_ids, "Age_Acceleration"])
AA_match$participant_id <- rownames(AA_match)

heat_dat <- clinical_hiv_prepared %>%
  select(participant_id, all_of(predictors)) %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  inner_join(AA_match, by = "participant_id")
y_vars <- c(names(AA_filter), "Multi_omics_AA")

# ---- Univariate regressions: age acceleration ~ clinical variable ----
get_category <- function(pred, level = NA_character_) {
  if (is.na(level)) return(pred)
  if (pred == "smoking_status") return(paste0("smoking_status_", level, "_vs_Never"))
  if (pred == "sex_at_birth") return(paste0("sex_at_birth_", level, "_vs_Male"))
  if (pred %in% c("reduced_egfr_lt90", "confirmed_diabetes", "established_ascvd",
                  "consistent_unquantifiable_vl_3y", "prior_mono_dual_art")) return(paste0(pred, "_", level, "_vs_0"))
  paste0(pred, "_", level)
}

heat_res_list <- list()
for (aa_col in y_vars) {
  for (pred in predictors) {
    d <- data.frame(y = heat_dat[[aa_col]], x = heat_dat[[pred]]) %>% filter(complete.cases(.))
    n_used <- nrow(d)
    if (n_used < 30 || length(unique(d$x)) < 2) next
    if (is.character(d$x) || is.factor(d$x)) d$x <- droplevels(factor(d$x))
    if (pred == "smoking_status" && "Never" %in% levels(d$x)) d$x <- relevel(d$x, ref = "Never")
    if (pred == "sex_at_birth" && "Male" %in% levels(d$x)) d$x <- relevel(d$x, ref = "Male")

    fit <- tryCatch(lm(y ~ x, data = d), error = function(e) NULL)
    if (is.null(fit)) next
    cf <- summary(fit)$coefficients
    ci <- confint.default(fit)

    if (is.factor(d$x)) {
      p_global <- tryCatch(drop1(fit, test = "F")["x", "Pr(>F)"], error = function(e) NA_real_)
      for (term in setdiff(rownames(cf), "(Intercept)")) {
        level <- sub("^x", "", term)
        heat_res_list[[length(heat_res_list) + 1]] <- data.frame(
          AA = aa_col, Feature = pred, Term = level, Category = get_category(pred, level),
          estimate = cf[term, "Estimate"], SE = cf[term, "Std. Error"],
          CI_low = ci[term, 1], CI_high = ci[term, 2],
          P_term = cf[term, "Pr(>|t|)"], P_global = p_global, N = n_used, stringsAsFactors = FALSE)
      }
    } else {
      heat_res_list[[length(heat_res_list) + 1]] <- data.frame(
        AA = aa_col, Feature = pred, Term = pred, Category = pred,
        estimate = cf["x", "Estimate"], SE = cf["x", "Std. Error"],
        CI_low = ci["x", 1], CI_high = ci["x", 2],
        P_term = cf["x", "Pr(>|t|)"], P_global = cf["x", "Pr(>|t|)"], N = n_used, stringsAsFactors = FALSE)
    }
  }
}

# FDR over all tests (SBP, HIV and ART duration not shown)
heat_res <- bind_rows(heat_res_list) %>%
  filter(!Feature %in% c("mean_sbp_mmHg", "hiv_duration_years", "art_duration_years")) %>%
  mutate(FDR_all = p.adjust(P_term, method = "BH"), FDR_plot = FDR_all)

# ---- Heatmap ----
fdr_levels <- c("not significant", "<0.05 negative", "<0.01 negative", "<0.001 negative", "<0.0001 negative",
                "<0.05 positive", "<0.01 positive", "<0.001 positive", "<0.0001 positive")
fdr_cols <- c("not significant" = "#FFFFFF",
              "<0.05 negative" = "#C6DBEF", "<0.01 negative" = "#6BAED6", "<0.001 negative" = "#3182BD", "<0.0001 negative" = "#08519C",
              "<0.05 positive" = "#FCBBA1", "<0.01 positive" = "#FC9272", "<0.001 positive" = "#DE2D26", "<0.0001 positive" = "#99000D")

heat_res <- heat_res %>%
  mutate(fdr_bin = factor(case_when(
    FDR_plot >= 0.05 ~ "not significant",
    estimate < 0 & FDR_plot < 0.0001 ~ "<0.0001 negative",
    estimate < 0 & FDR_plot < 0.001  ~ "<0.001 negative",
    estimate < 0 & FDR_plot < 0.01   ~ "<0.01 negative",
    estimate < 0                     ~ "<0.05 negative",
    FDR_plot < 0.0001 ~ "<0.0001 positive",
    FDR_plot < 0.001  ~ "<0.001 positive",
    FDR_plot < 0.01   ~ "<0.01 positive",
    TRUE              ~ "<0.05 positive"), levels = fdr_levels))

category_levels <- c("sex_at_birth_Female_vs_Male", "smoking_status_Former_vs_Never", "smoking_status_Current_vs_Never",
                     "BMI", "mean_dbp_mmHg", "reduced_egfr_lt90_1_vs_0", "fibroscan_cap_median",
                     "confirmed_diabetes_1_vs_0", "established_ascvd_1_vs_0", "cd4_nadir_per100", "cd4_latest_per100",
                     "cd8_latest_per100", "CD4_CD8_ratio", "consistent_unquantifiable_vl_3y_1_vs_0", "prior_mono_dual_art_1_vs_0")
category_levels <- c(category_levels[category_levels %in% heat_res$Category],
                     setdiff(unique(heat_res$Category), category_levels))
x_lab <- c(sex_at_birth_Female_vs_Male = "Sex", smoking_status_Former_vs_Never = "Former smoking",
           smoking_status_Current_vs_Never = "Current smoking", BMI = "BMI", mean_dbp_mmHg = "DBP",
           reduced_egfr_lt90_1_vs_0 = "eGFR <90", fibroscan_cap_median = "FibroScan CAP",
           confirmed_diabetes_1_vs_0 = "Diabetes", established_ascvd_1_vs_0 = "ASCVD",
           cd4_nadir_per100 = "CD4 nadir (per 100 cells/µL)", cd4_latest_per100 = "Latest CD4 (per 100 cells/µL)",
           cd8_latest_per100 = "Latest CD8 (per 100 cells/µL)", CD4_CD8_ratio = "CD4/CD8",
           consistent_unquantifiable_vl_3y_1_vs_0 = "Sustained VL suppression", prior_mono_dual_art_1_vs_0 = "Prior mono/dual ART")
aa_lab <- c(Multi_omics_AA = "Multi-omics AA", EAA_Grimage = "GrimAge EAA", EAA_Phenoage = "PhenoAge EAA",
            EAA_Hannum = "Hannum EAA", EAA_Horvath = "Horvath EAA")

heat_res <- heat_res %>% mutate(Category = factor(Category, levels = category_levels), AA = factor(AA, levels = y_vars))

p_heat <- ggplot(heat_res, aes(Category, AA, fill = fdr_bin)) +
  geom_tile(colour = "grey85", linewidth = 0.2) +
  scale_fill_manual(values = fdr_cols, limits = fdr_levels, drop = FALSE, name = "FDR") +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(x_lab), x_lab[x], x)) +
  scale_y_discrete(labels = function(x) ifelse(x %in% names(aa_lab), aa_lab[x], x)) +
  theme_minimal(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.title = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"))

write.csv(heat_res, file.path(out_dir, "tabpfn_heatmap_AA_multiomics_clinical_FDR3.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "tabpfn_heatmap_AA_multiomics_clinical_FDR3.pdf"), p_heat,
       width = max(170, 8 * length(category_levels)), height = max(55, 8 * length(y_vars)), units = "mm")
