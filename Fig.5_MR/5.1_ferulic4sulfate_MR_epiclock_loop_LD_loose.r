library(data.table)
library(dplyr)
library(remotes)
library(TwoSampleMR)
library(ggplot2)
library(stringr)
library(ieugwasr)

packageVersion("TwoSampleMR")
R.version.string

# load("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD/MR_ferulic_4_sulfate.Rdata")

load("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD/MR_ferulic_4_sulfate_correctLD.Rdata")

outcomes_list <- list(
  list(
    name = "PhenoAge_AA",
    file = "/vol/projects/yzhang/500FG_aging/output/16_MR/AA_qtl_filter.RDS",
    sample_size = 34710
  ),
  list(
    name = "Hannum_AA",
    file = "/vol/projects/yzhang/500FG_aging/input/MR/Hannum_qtl_filter.RDS",  
    sample_size = 34710 
  ),
  list(
    name = "IEAA",
    file = "/vol/projects/yzhang/500FG_aging/input/MR/IEAA_qtl_filter.RDS",    
    sample_size = 34710 
  ),
  list(
    name = "GrimAge_AA",
    file = "/vol/projects/yzhang/500FG_aging/input/MR/GrimAge_qtl_filter.RDS", 
    sample_size = 34710
  )
)

N_EXP <- 2392   # exposure sample size
F_MIN <- 10      # minimum F-stat threshold


# Check optional packages
has_mrpresso <- requireNamespace("MRPRESSO", quietly = TRUE)
has_raps     <- requireNamespace("mr.raps", quietly = TRUE)
has_grip     <- any("mr_grip" %in% mr_method_list()$obj)

# Helper functions
safe1 <- function(x, col) {
  if (is.null(x) || !nrow(x) || !(col %in% names(x))) return(NA_real_)
  as.numeric(x[[col]][1])
}

pickrow <- function(df, pat) {
  if (is.null(df) || !nrow(df)) return(NULL)
  df[grepl(pat, df$method), , drop = FALSE] %>% head(1)
}


NHW_mqtl <- fread("/vol/projects/yzhang/500FG_aging/input/MR/ferulic_4_sulfate/phenocode-X100005389_Plasma_NHW.tsv.gz")  ##2392 samples #download from 
head(NHW_mqtl)
dim(NHW_mqtl)

rm(NHW_mqtl)

NHW_mqtl_sig <- NHW_mqtl %>%
  filter(pval < 1e-06)
head(NHW_mqtl_sig)
dim(NHW_mqtl_sig)

write.csv(NHW_mqtl_sig,"/vol/projects/yzhang/500FG_aging/output/16_MR/NHW_mqtl_sig1e6.csv",row.names=TRUE)

##check the LD pruning 
snps_ind1 <- read.table("/vol/projects/mballan/send/EUR_pruned.prune.in",  #snps ID
                        header = FALSE, 
                        stringsAsFactors = FALSE)
dim(snps_ind1)
head(snps_ind1)
# Overlap数量
length(intersect(snps_ind1$V1, NHW_mqtl$rsids))

# 只在snps_ind1中有的
length(setdiff(snps_ind1$V1, NHW_mqtl$rsids))

# 只在NHW_mqtl中有的
length(setdiff(NHW_mqtl$rsids, snps_ind1$V1))

###using TwoSampleMR for ld pruning
# 2. Format as exposure data
exp_dat <- NHW_mqtl_sig %>%
  as.data.frame() %>%
  format_data(
    type = "exposure",
    snp_col = "rsids", 
    beta_col = "beta", 
    se_col = "sebeta", 
    effect_allele_col = "ref",      
    other_allele_col = "alt",       
    pval_col = "pval",
    chr_col = "chrom", 
    pos_col = "pos"
  ) %>%
  mutate(exposure = "Ferulic 4 sulfate")

if (!"eaf.exposure" %in% names(exp_dat)) exp_dat$eaf.exposure <- NA_real_

message("Exposure data formatted: ", nrow(exp_dat), " SNPs")

# ============================================================
# 2. Load 500FG reference and extract rsIDs
# ============================================================
message("Loading 500FG reference panel...")

bim <- fread("/vol/projects/CIIM/meta_cQTL/out/500FG/genotype/allchr.bim", 
             header = FALSE,
             col.names = c("chr", "snp_id", "cm", "pos", "a1", "a2"))

# Extract rsID from the complex ID format (after semicolon)
# Format: chr1:794707:T:C;rs148120343 -> rs148120343
bim <- bim %>%
  mutate(rsid = ifelse(grepl(";", snp_id), 
                       sub(".*;", "", snp_id),  # Extract part after semicolon
                       snp_id))

message("Reference panel loaded: ", nrow(bim), " SNPs")

# ============================================================
# 3. Check overlap with extracted rsIDs
# ============================================================
n_overlap <- sum(exp_dat$SNP %in% bim$rsid)
message("Overlapping SNPs after rsID extraction: ", n_overlap, " / ", nrow(exp_dat))

# ============================================================
# 4. Create temporary .bim file with rsIDs for PLINK
# ============================================================
message("Creating temporary reference panel with rsIDs...")

# Create temp directory
temp_dir <- tempdir()
temp_prefix <- file.path(temp_dir, "ref_rsid")

# Write new .bim file with rsIDs
bim_new <- bim %>%
  mutate(V2 = rsid) %>%
  select(chr, V2, cm, pos, a1, a2)

fwrite(bim_new, paste0(temp_prefix, ".bim"), 
       sep = "\t", col.names = FALSE)

# Create symlinks for .bed and .fam files
file.symlink("/vol/projects/CIIM/meta_cQTL/out/500FG/genotype/allchr.bed",
             paste0(temp_prefix, ".bed"))
file.symlink("/vol/projects/CIIM/meta_cQTL/out/500FG/genotype/allchr.fam",
             paste0(temp_prefix, ".fam"))

message("Temporary reference created: ", temp_prefix)

# ============================================================
# 5. LD Clumping with corrected reference
# ============================================================
message("Performing LD clumping...")

library(ieugwasr)

clump_input <- data.frame(
  rsid = exp_dat$SNP,
  pval = exp_dat$pval.exposure
)

clumped_snps <- ld_clump(
  dat = clump_input,
  clump_kb = 5000,
  clump_r2 = 0.5,
  clump_p = 1e-05,
  bfile = temp_prefix,
  plink_bin = "/vol/projects/yzhang/500FG_aging/code/prediction/prediction_prs/plink"
)

# ============================================================
# 6. Filter to keep independent SNPs
# ============================================================
exp_dat_clumped <- exp_dat %>%
  filter(SNP %in% clumped_snps$rsid)

message("Independent SNPs after LD clumping: ", nrow(exp_dat_clumped))

# View results
print(exp_dat_clumped)

# Load independent SNPs  ！！！！！may not correct!!!!!
# snps_ind <- read.table("/vol/projects/mballan/send/EUR_pruned.prune.in",  #snps ID
#                        header = FALSE, 
#                        stringsAsFactors = FALSE)
# snps_ind <- read.table("/vol/projects/mballan/send/pos.txt",  #snps postion
#                        header = FALSE, 
#                        stringsAsFactors = FALSE)
# snps_ind <- read.table("/vol/projects/mballan/send/pruned.03.prune.in",  #snps 0.3 treshold
#                        header = FALSE, 
#                        stringsAsFactors = FALSE)
# head(snps_ind)
# # got position info
# NHW_mqtl_sig <- NHW_mqtl_sig %>%
#   mutate(pos_id = paste0(chrom, ":", pos))

# # Filter for independent SNPs
# ferulic_sig_common <- NHW_mqtl_sig %>%
#   filter(rsids %in% snps_ind$V1 | pos_id %in% snps_ind$V1)
# head(ferulic_sig_common)

# message("  ✓ Total exposure SNPs (p<5e-08): ", nrow(NHW_mqtl_sig))
# message("  ✓ Independent SNPs after LD pruning: ", nrow(ferulic_sig_common))
# message("    - Matched by rsID: ", sum(ferulic_sig_common$rsids %in% snps_ind$V1))
# message("    - Matched by position: ", sum(ferulic_sig_common$pos_id %in% snps_ind$V1))

# # Format exposure data 
# exp_dat <- ferulic_sig_common %>%
#   as.data.frame() %>%
#   format_data(
#     type = "exposure",
#     snp_col = "rsids", 
#     beta_col = "beta", 
#     se_col = "sebeta", 
#     effect_allele_col = "ref",      
#     other_allele_col = "alt",       
#     pval_col = "pval",
#     chr_col = "chrom", 
#     pos_col = "pos"
#   ) %>%
#   mutate(exposure = "Feruli 4 sulfate")

# if (!"eaf.exposure" %in% names(exp_dat)) exp_dat$eaf.exposure <- NA_real_

# message("  ✓ Exposure data formatted: ", nrow(exp_dat), " SNPs")

# LOOP THROUGH ALL OUTCOMES
# ============================================================
# set output path
output_base_dir <- "/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD_loose"

# Initialize results storage
all_results <- list()
all_mr_res <- list()

for (i in seq_along(outcomes_list)) {
  
  outcome_info <- outcomes_list[[i]]
  outcome_name <- outcome_info$name
  outcome_file <- outcome_info$file
  N_OUT <- outcome_info$sample_size
  
  
  message("PROCESSING OUTCOME ", i, "/", length(outcomes_list), ": ", outcome_name)
 
  
  # -------------------- Load Outcome Data --------------------
  tryCatch({
    
    message("Loading outcome data: ", outcome_file)
    outcome_data <- readRDS(outcome_file)
    message("  ✓ Loaded ", nrow(outcome_data), " rows")
    
    # Format outcome data
    out_dat <- outcome_data %>%
      rename(
        effect_allele = A1, 
        other_allele = A2, 
        beta = Effect, 
        se = SE, 
        pval = P
      ) %>%
      format_data(
        type = "outcome",
        snp_col = "SNP", 
        beta_col = "beta",
        se_col = "se", 
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele", 
        pval_col = "pval",
        chr_col = "chr", 
        pos_col = "pos"
      ) %>%
      mutate(outcome = outcome_name)
    
    if (!"eaf.outcome" %in% names(out_dat)) out_dat$eaf.outcome <- NA_real_
    
    message("  ✓ Outcome SNPs: ", nrow(out_dat))
    
    # -------------------- Harmonise --------------------
    message("Harmonising...")
    
    dat_all <- harmonise_data(exp_dat_clumped, out_dat, action = 2) ####use clump data
    
    # Robustly coerce list columns
    if ("remove" %in% names(dat_all) && is.list(dat_all$remove)) {
      dat_all$remove <- unlist(dat_all$remove)
    }
    if ("ambiguous" %in% names(dat_all) && is.list(dat_all$ambiguous)) {
      dat_all$ambiguous <- unlist(dat_all$ambiguous)
    }
    
    dat_all <- dat_all %>%
      filter(remove == FALSE, is.na(ambiguous) | ambiguous == FALSE)
    
    if (nrow(dat_all) == 0) {
      message("   No overlapping SNPs after harmonisation. Skipping ", outcome_name)
      next
    }
    
    message("  ✓ Harmonised SNPs: ", nrow(dat_all))
    
    dat_all$samplesize.exposure <- N_EXP
    dat_all$samplesize.outcome  <- N_OUT
    
    # -------------------- Calculate F-statistics --------------------
    message("Calculating F-statistics...")
    
    dat_all <- dat_all %>% 
      mutate(F_stat = (beta.exposure^2)/(se.exposure^2))
    
    F_mean <- mean(dat_all$F_stat, na.rm = TRUE)
    message("  ✓ Mean F-statistic: ", round(F_mean, 2))
    
    dat_filtered <- if (F_MIN > 0) {
      before <- nrow(dat_all)
      dat_f <- dat_all %>% filter(F_stat >= F_MIN)
      after <- nrow(dat_f)
      message("  ✓ Filtered SNPs with F >= ", F_MIN, ": ", after, " (removed ", before - after, ")")
      dat_f
    } else {
      dat_all
    }
    
    if (nrow(dat_filtered) == 0) {
      message("  No SNPs remaining after F-statistic filtering. Skipping ", outcome_name)
      next
    }
    
    # -------------------- Run MR Analyses --------------------
    message("Running MR analyses...")
    
    ns <- length(unique(dat_filtered$SNP))
    message("  ✓ Number of SNPs used: ", ns)
    
    # Define methods list
    mlist <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median")
    if (has_raps) mlist <- c(mlist, "mr_raps")
    if (has_grip) mlist <- c(mlist, "mr_grip")
    
    # Run MR
    mr_res <- tryCatch(mr(dat_filtered, method_list = mlist), error = function(e) {
      message("  MR failed: ", e$message)
      return(NULL)
    })
    
    if (is.null(mr_res)) {
      message("  MR analysis failed. Skipping ", outcome_name)
      next
    }
    
    # Extract primary result
    wald <- pickrow(mr_res, "^Wald ratio$")
    ivw  <- pickrow(mr_res, "Inverse variance weighted")
    
    if (ns == 1 && !is.null(wald)) {
      method_main <- "Wald ratio"
      beta_main <- safe1(wald, "b")
      se_main   <- safe1(wald, "se")
      p_main    <- safe1(wald, "pval")
    } else {
      method_main <- "Inverse variance weighted (random effects)"
      beta_main <- safe1(ivw, "b")
      se_main   <- safe1(ivw, "se")
      p_main    <- safe1(ivw, "pval")
    }
    
    message("  ✓ Primary MR: Beta = ", round(beta_main, 4), 
            ", SE = ", round(se_main, 4), 
            ", P = ", format(p_main, scientific = TRUE, digits = 3))
    
    # -------------------- Heterogeneity Test --------------------
    message("Testing heterogeneity...")
    
    het <- tryCatch(mr_heterogeneity(dat_filtered), error = function(e) NULL)
    Q_p_ivw <- NA_real_
    
    if (!is.null(het) && nrow(het)) {
      hi <- het[het$method == "Inverse variance weighted", , drop = FALSE]
      if (nrow(hi)) {
        Q_p_ivw <- as.numeric(hi$Q_pval[1])
        message("  ✓ Cochran's Q p-value: ", format(Q_p_ivw, scientific = TRUE, digits = 3))
      }
    }
    
    # -------------------- Pleiotropy Test --------------------
    message(" Testing pleiotropy...")
    
    pleio <- tryCatch(mr_pleiotropy_test(dat_filtered), error = function(e) NULL)
    egger_intercept_p <- if (!is.null(pleio) && nrow(pleio)) as.numeric(pleio$pval[1]) else NA_real_
    
    if (!is.na(egger_intercept_p)) {
      message("  ✓ MR-Egger p-value: ", format(egger_intercept_p, scientific = TRUE, digits = 3))
    }
    
    # -------------------- MR-PRESSO --------------------
    message(" Running MR-PRESSO...")
    
    presso_p <- NA_real_
    presso_results <- NULL
    
    if (has_mrpresso && ns >= 3) {
      pr <- tryCatch(
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", 
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",   
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE,
          data = as.data.frame(dat_filtered), 
          NbDistribution = 1000, 
          SignifThreshold = 0.05
        ),
        error = function(e) NULL
      )
      
      if (!is.null(pr) && !is.null(pr$`MR-PRESSO results`)) {
        presso_results <- pr
        presso_p <- tryCatch(pr$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
        if (!is.na(presso_p)) {
          message("  ✓ MR-PRESSO p-value: ", format(presso_p, scientific = TRUE, digits = 3))
        }
      }
    } else {
      message("   MR-PRESSO: Need ≥3 SNPs (have ", ns, ")")
    }
    
    # -------------------- Extract Optional Methods --------------------
    raps <- pickrow(mr_res, "RAPS")
    grip <- pickrow(mr_res, "GRIP")
    
    # -------------------- Compile Results --------------------
    results_table <- tibble::tibble(
      exposure = "Feruli 4 sulfate",
      outcome = outcome_name,
      n_SNPs_used = ns,
      F_mean_exposure_strength = F_mean,
      Q_pval_cochrans_Q_IVW = Q_p_ivw,
      Egger_intercept_pval_pleiotropy = egger_intercept_p,
      MR_PRESSO_global_pval_outlier_test = presso_p,
      MR_method_primary = method_main,
      MR_beta_estimate = beta_main,
      MR_standard_error = se_main,
      MR_pval_main = p_main,
      MR_RAPS_beta = safe1(raps, "b"),
      MR_RAPS_pval = safe1(raps, "pval"),
      MR_GRIP_beta = safe1(grip, "b"),
      MR_GRIP_pval = safe1(grip, "pval")
    )
    
    # Store results
    all_results[[outcome_name]] <- results_table
    all_mr_res[[outcome_name]] <- list(
      mr_res = mr_res,
      het = het,
      pleio = pleio,
      presso = presso_results,
      dat_filtered = dat_filtered
    )
    
    # -------------------- Save Individual Results --------------------
    output_dir <- file.path(output_base_dir, paste0("MR_output_", outcome_name))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save results
    data.table::fwrite(results_table, 
                       file.path(output_dir, "MR_sensitivity_results.csv"), 
                       sep = ";", dec = ",", quote = FALSE, na = "NA")
    
    if (!is.null(mr_res)) {
      write.csv(mr_res, file.path(output_dir, "MR_all_methods.csv"), row.names = FALSE)
    }
    if (!is.null(het)) {
      write.csv(het, file.path(output_dir, "MR_heterogeneity_detailed.csv"), row.names = FALSE)
    }
    if (!is.null(pleio)) {
      write.csv(pleio, file.path(output_dir, "MR_pleiotropy_detailed.csv"), row.names = FALSE)
    }
    if (!is.null(presso_results)) {
      saveRDS(presso_results, file.path(output_dir, "MR_PRESSO_full_results.rds"))
    }
    
    write.csv(dat_filtered, file.path(output_dir, "harmonised_data_with_F.csv"), row.names = FALSE)
    
    # -------------------- Generate Plots --------------------
    message(" Generating plots...")
    
    plot_dir <- file.path(output_base_dir, paste0("MR_plots_", outcome_name))
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Scatter plot
    tryCatch({
      p1 <- mr_scatter_plot(mr_res, dat_filtered)
      ggsave(file.path(plot_dir, "scatter_plot.png"), p1[[1]], 
             width = 10, height = 8, dpi = 300)
    }, error = function(e) message("  ⚠️ Scatter plot failed"))
    
    # Forest plot
    tryCatch({
      single_snp_res <- mr_singlesnp(dat_filtered)
      p2 <- mr_forest_plot(single_snp_res)
      ggsave(file.path(plot_dir, "forest_plot.png"), p2[[1]], 
             width = 10, height = 12, dpi = 300)
    }, error = function(e) message("  ⚠️ Forest plot failed"))
    
    # Leave-one-out plot
    if (ns >= 3) {
      tryCatch({
        loo_res <- mr_leaveoneout(dat_filtered)
        p3 <- mr_leaveoneout_plot(loo_res)
        ggsave(file.path(plot_dir, "leaveoneout_plot.png"), p3[[1]], 
               width = 10, height = 12, dpi = 300)
      }, error = function(e) message("  ⚠️ Leave-one-out plot failed"))
    }
    
    # Funnel plot
    if (ns >= 3) {
      tryCatch({
        single_snp_res <- mr_singlesnp(dat_filtered)
        p4 <- mr_funnel_plot(single_snp_res)
        ggsave(file.path(plot_dir, "funnel_plot.png"), p4[[1]], 
               width = 8, height = 8, dpi = 300)
      }, error = function(e) message("  ⚠️ Funnel plot failed"))
    }
    
    message("  Plots saved in: ", plot_dir)
    
    # -------------------- Display Summary --------------------
    cat("\n")
    cat("="*60, "\n")
    cat(" SUMMARY: ", outcome_name, "\n")
    cat("="*60, "\n")
    cat("SNPs used:            ", ns, "\n")
    cat("Mean F-statistic:     ", round(F_mean, 2), "\n")
    cat("Beta (IVW):           ", round(beta_main, 4), "\n")
    cat("P-value:              ", format(p_main, scientific = TRUE, digits = 3), "\n")
    cat("Heterogeneity (Q):    ", if (is.na(Q_p_ivw)) "NA" else format(Q_p_ivw, scientific = TRUE, digits = 3), "\n")
    cat("Pleiotropy (Egger):   ", if (is.na(egger_intercept_p)) "NA" else format(egger_intercept_p, scientific = TRUE, digits = 3), "\n")
    cat("Outliers (PRESSO):    ", if (is.na(presso_p)) "NA" else format(presso_p, scientific = TRUE, digits = 3), "\n")
    cat("\n")
    
    message(" Completed: ", outcome_name)
    
  }, error = function(e) {
    message(" ERROR processing ", outcome_name, ": ", e$message)
  })
  
}

# ============================================================
#  COMBINE ALL RESULTS
# ============================================================


message(" COMBINING ALL RESULTS")


if (length(all_results) > 0) {
  
  combined_results <- bind_rows(all_results)
  
  # Save combined results
  combined_output_dir <- file.path(output_base_dir, "MR_combined_results")
  dir.create(combined_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(combined_results, 
            file.path(combined_output_dir, "MR_all_outcomes_combined.csv"), 
            row.names = FALSE)
  
  data.table::fwrite(combined_results, 
                     file.path(combined_output_dir, "MR_all_outcomes_combined_EUR.csv"), 
                     sep = ";", dec = ",", quote = FALSE, na = "NA")
  
  message(" Combined results saved: ", file.path(combined_output_dir, "MR_all_outcomes_combined.csv"))
  
  # Display combined summary
  cat("\n")
  cat(" COMBINED RESULTS SUMMARY\n")
  
  print(combined_results %>% 
    select(outcome, n_SNPs_used, MR_beta_estimate, MR_standard_error, MR_pval_main) %>%
    mutate(
      CI_lower = MR_beta_estimate - 1.96 * MR_standard_error,
      CI_upper = MR_beta_estimate + 1.96 * MR_standard_error
    ))
  
} else {
  message(" No results to combine!")
}


message(" ALL ANALYSES COMPLETE!")
message("\nProcessed ", length(all_results), " out of ", length(outcomes_list), " outcomes")
message("\nResults saved in: ", output_base_dir)

#  COMBINE ALL RESULTS
# Combined Forest Plot
  tryCatch({
    forest_data <- combined_results %>%
      mutate(
        CI_lower = MR_beta_estimate - 1.96 * MR_standard_error,
        CI_upper = MR_beta_estimate + 1.96 * MR_standard_error
      )
    
    p_forest <- ggplot(forest_data, aes(x = MR_beta_estimate, y = outcome)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
      geom_point(size = 3) +
      labs(x = "Beta Estimate (95% CI)", y = "Outcome") +
      theme_minimal()
    
   ggsave(file.path(combined_output_dir, "combined_forest_plot.png"), 
          p_forest, width = 10, height = 6, dpi = 300)
    
    message(" Forest plot saved")
  }, error = function(e) message(" Forest plot failed: ", e$message))
  
  # Combined Scatter Plot
  tryCatch({
    scatter_data <- bind_rows(lapply(names(all_mr_res), function(name) {
      cbind(all_mr_res[[name]]$dat_filtered, outcome = name)
    }))
    
    p_scatter <- ggplot(scatter_data, aes(x = beta.exposure, y = beta.outcome)) +
      geom_point(aes(color = outcome)) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
      labs(x = "SNP effect on Exposure", y = "SNP effect on Outcome") +
      theme_bw() +
      theme(legend.position = "none")
    
   ggsave(file.path(combined_output_dir, "combined_scatter_plot.png"), 
          p_scatter, width = 12, height = 10, dpi = 300)
    
    message(" Scatter plot saved")
  }, error = function(e) message(" Scatter plot failed: ", e$message))
  

message(" ALL ANALYSES COMPLETE!")
message("\nProcessed ", length(all_results), " out of ", length(outcomes_list), " outcomes")
message("\nResults saved in: ", output_base_dir)

 # Combined Forest Plot with individual SNPs
  tryCatch({
    # Prepare data for all SNPs from all outcomes
    all_snp_data <- bind_rows(lapply(names(all_mr_res), function(outcome_name) {
      dat <- all_mr_res[[outcome_name]]$dat_filtered
      snp_effects <- dat %>%
        mutate(
          # Calculate MR effect size for each SNP (beta.outcome / beta.exposure)
          mr_effect = beta.outcome / beta.exposure,
          # Calculate SE using delta method
          se_mr = abs(mr_effect) * sqrt((se.outcome/beta.outcome)^2 + (se.exposure/beta.exposure)^2),
          CI_lower = mr_effect - 1.96 * se_mr,
          CI_upper = mr_effect + 1.96 * se_mr,
          label = SNP,
          outcome = outcome_name,
          type = "SNP"
        ) %>%
        select(outcome, label, beta = mr_effect, CI_lower, CI_upper, type)
      return(snp_effects)
    }))
    
    # Add IVW summary results with p-values
    summary_data <- combined_results %>%
      mutate(
        CI_lower = MR_beta_estimate - 1.96 * MR_standard_error,
        CI_upper = MR_beta_estimate + 1.96 * MR_standard_error,
        label = paste0("IVW (n=", n_SNPs_used, ")"),
        beta = MR_beta_estimate,
        pval = MR_pval_main,
        type = "Summary"
      ) %>%
      select(outcome, label, beta, CI_lower, CI_upper, pval, type)
      # Combine SNP and summary data
    all_snp_data$pval <- NA_real_  # Add pval column to SNP data
    
    # Define outcome order
    outcome_order <- c("PhenoAge_AA", "Hannum_AA", "GrimAge_AA", "IEAA")
    
    forest_data <- bind_rows(all_snp_data, summary_data) %>%
      mutate(outcome_group = factor(outcome, levels = outcome_order)) %>%
      arrange(outcome_group, type) %>%
      group_by(outcome_group) %>%
      mutate(y_position_in_group = row_number()) %>%
      ungroup() %>%
      mutate(
        # Add extra spacing between outcomes (5 units), tight within outcomes (0.8 units)
        outcome_number = as.numeric(outcome_group),
        y_position = y_position_in_group * 0.8 + (outcome_number - 1) * 5
      )
    
    # Create text labels for IVW results
    ivw_labels <- forest_data %>%
      filter(type == "Summary") %>%
      mutate(
        text_label = paste0("Beta=", round(beta, 3), 
                           ", p=", format(round(pval, 4), nsmall = 4))
      )
    
    # Calculate fixed x position for text (right side of plot)
    x_max <- max(forest_data$CI_upper, na.rm = TRUE)
    x_min <- min(forest_data$CI_lower, na.rm = TRUE)
    x_range <- x_max - x_min
    x_text_pos <- x_max + x_range * 0.05  # Text position (much closer)
    x_panel_max <- x_text_pos - x_range * 0.02  # Panel ends just before text
    x_limit <- x_max + x_range * 0.2     # Extend x-axis limit for text display
    
    # Create white rectangle data to hide grid lines on the right
    rect_data <- data.frame(
      xmin = x_panel_max,
      xmax = x_limit,
      ymin = -Inf,
      ymax = Inf,
      outcome_group = unique(forest_data$outcome_group)
    )
    
    # Create plot
    p_forest <- ggplot(forest_data, aes(x = beta, y = y_position)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = type),
                     height = 0.3, size = 0.7) +
      geom_point(aes(color = type, size = type, shape = type)) +
      geom_rect(data = rect_data, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "white", color = NA, inherit.aes = FALSE) +
      geom_label(data = ivw_labels, 
                 aes(x = x_text_pos, y = y_position, label = text_label),
                 hjust = 0, size = 3.5, fontface = "bold", 
                 color = "black", fill = "white", 
                 label.padding = unit(0.25, "lines"),
                 label.size = NA) +
      scale_color_manual(values = c("SNP" = "steelblue", "Summary" = "#E63946")) +
      scale_size_manual(values = c("SNP" = 2, "Summary" = 4)) +
      scale_shape_manual(values = c("SNP" = 16, "Summary" = 18)) +
      scale_x_continuous(limits = c(x_min, x_limit), expand = c(0.01, 0)) +
      scale_y_continuous(breaks = forest_data$y_position,
                        labels = forest_data$label) +
      coord_cartesian(xlim = c(x_min, x_limit), clip = "off") +
      facet_grid(outcome_group ~ ., scales = "free_y", space = "free_y") +
      labs(title = "MR Forest Plot: Ferulic 4 sulfate  -> Aging Markers",
           x = "MR effect size", 
           y = "") +
      theme_minimal(base_size = 11) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        axis.text.y = element_text(size = 9),
        strip.text.y = element_text(face = "bold", size = 10, margin = margin(0, 2, 0, 2)),
        strip.background = element_rect(fill = "lightblue", color = "gray40"),
        strip.placement = "outside",
        legend.position = "bottom",
        panel.spacing.y = unit(1, "lines"),
        plot.margin = margin(5, 25, 5, 5),
        axis.title.x = element_text(margin = margin(t = 5))
      )
    
    ggsave(file.path(combined_output_dir, "combined_forest_plot_detailed_LD.pdf"), 
           p_forest, width = 13, height = 8, bg = "white")
    
    message("Detailed forest plot saved (with all SNPs)")
      }, error = function(e) message("Forest plot failed: ", e$message))

save.image("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD/MR_ferulic_4_sulfate.Rdata")

save.image("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD/MR_ferulic_4_sulfate_correctLD.Rdata")

load("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate_LD/MR_ferulic_4_sulfate_correctLD.Rdata")

##nature sytle plot
# Combined Forest Plot — Nature style (85 mm × 60 mm, font 5 pt)
tryCatch({
  # Prepare data for all SNPs from all outcomes
  all_snp_data <- bind_rows(lapply(names(all_mr_res), function(outcome_name) {
    dat <- all_mr_res[[outcome_name]]$dat_filtered
    snp_effects <- dat %>%
      mutate(
        mr_effect = beta.outcome / beta.exposure,
        se_mr = abs(mr_effect) * sqrt((se.outcome/beta.outcome)^2 + (se.exposure/beta.exposure)^2),
        CI_lower = mr_effect - 1.96 * se_mr,
        CI_upper = mr_effect + 1.96 * se_mr,
        label = SNP,
        outcome = outcome_name,
        type = "SNP"
      ) %>%
      select(outcome, label, beta = mr_effect, CI_lower, CI_upper, type)
    return(snp_effects)
  }))

  summary_data <- combined_results %>%
    mutate(
      CI_lower = MR_beta_estimate - 1.96 * MR_standard_error,
      CI_upper = MR_beta_estimate + 1.96 * MR_standard_error,
      label = paste0("IVW (n=", n_SNPs_used, ")"),
      beta = MR_beta_estimate,
      pval = MR_pval_main,
      type = "Summary"
    ) %>%
    select(outcome, label, beta, CI_lower, CI_upper, pval, type)

  all_snp_data$pval <- NA_real_

  outcome_order <- c("PhenoAge_AA", "Hannum_AA", "GrimAge_AA", "IEAA")
  outcome_labels <- c("PhenoAge_AA" = "PhenoAge AA", "Hannum_AA" = "Hannum AA",
                     "GrimAge_AA" = "GrimAge AA", "IEAA" = "IEAA")
    
  forest_data <- bind_rows(all_snp_data, summary_data) %>%
    mutate(outcome_group = factor(outcome, levels = outcome_order)) %>%
    arrange(outcome_group, type) %>%
    group_by(outcome_group) %>%
    mutate(y_position_in_group = row_number()) %>%
    ungroup() %>%
    mutate(
      outcome_number = as.numeric(outcome_group),
      y_position = y_position_in_group * 0.6 + (outcome_number - 1) * 4
    )

  ivw_labels <- forest_data %>%
    filter(type == "Summary") %>%
    mutate(
          text_label = paste0("Beta=", round(beta, 3), "\n", "p=", format(round(pval, 3), nsmall = 3))
    )

  x_max <- max(forest_data$CI_upper, na.rm = TRUE)
  x_min <- min(forest_data$CI_lower, na.rm = TRUE)
  x_range <- x_max - x_min
  x_text_pos <- x_max + x_range * 0.08
  x_panel_max <- x_text_pos - x_range * 0.02
  x_limit <- x_max + x_range * 0.35

  rect_data <- data.frame(
    xmin = x_panel_max,
    xmax = x_limit,
    ymin = -Inf,
    ymax = Inf,
    outcome_group = unique(forest_data$outcome_group)
  )

  # Nature-style colours (colourblind-friendly, print-safe)
  col_snp     <- "#2166AC"
  col_summary <- "#B2182B"

  p_forest <- ggplot(forest_data, aes(x = beta, y = y_position)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.25) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = type),
                   height = 0.2, linewidth = 0.35) +
    geom_point(aes(color = type, size = type, shape = type)) +
    geom_rect(data = rect_data,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = NA, inherit.aes = FALSE) +
    geom_label(data = ivw_labels,
               aes(x = x_text_pos, y = y_position, label = text_label),
               hjust = 0, size = 1.35,
               color = "black", fill = "white",
               label.padding = unit(0.12, "lines"),
               label.size = NA, label.r = unit(0, "lines")) +
    scale_color_manual(values = c("SNP" = col_snp, "Summary" = col_summary)) +
    scale_size_manual(values = c("SNP" = 0.8, "Summary" = 1.4)) +
    scale_shape_manual(values = c("SNP" = 16, "Summary" = 18)) +
    scale_x_continuous(limits = c(x_min, x_limit), expand = c(0.02, 0)) +
    scale_y_continuous(breaks = forest_data$y_position, labels = forest_data$label) +
    coord_cartesian(xlim = c(x_min, x_limit), clip = "off") +
       facet_grid(outcome_group ~ ., scales = "free_y", space = "free_y",
           labeller = as_labeller(outcome_labels)) +
     labs(#title = "Ferulic 4 sulfate → Aging markers",
         x = "MR effect size",
         y = NULL) +
    theme_minimal(base_size = 5) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.2),
      axis.text = element_text(size = 5, colour = "black"),
      axis.title = element_text(size = 5, colour = "black"),
      axis.title.x = element_text(margin = margin(t = 1)),
      plot.title = element_text(size = 5, hjust = 0.5, margin = margin(b = 1)),
      strip.text.y = element_text( size = 5, colour = "black",
                                  margin = margin(0, 0.5, 0, 0.5)),
      strip.background = element_rect(fill = "grey96", color = NA),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.text = element_text(size = 5, colour = "black"),
      legend.title = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.spacing.x = unit(0.15, "cm"),
      panel.spacing.y = unit(0.4, "lines"),
      plot.margin = margin(4, 12, 2, 0),
      axis.ticks = element_line(linewidth = 0.2, colour = "black"),
      axis.ticks.length = unit(0.05, "cm")
    )

  # Nature single column: 85 mm × 60 mm
  fig_w <- 85 / 25.4   # inches
  fig_h <- 62 / 25.4

  ggsave(file.path(combined_output_dir, "combined_forest_plot_detailed_LD_nature.pdf"),
         p_forest, width = fig_w, height = fig_h, units = "in", dpi = 300,
          bg = "white"#, device = cairo_pdf
        )

  message("Detailed forest plot saved (Nature style, 85×60 mm)")
}, error = function(e) message("Forest plot failed: ", e$message))


# ============================================================
# REVERSE MR ANALYSIS: Epigenetic Clocks → Feruli 4 sulfate
# ============================================================
# This tests whether aging phenotypes causally affect metabolite levels
# (to rule out reverse causation)

message("=" , strrep("=", 60))
message(" STARTING REVERSE MR ANALYSIS")
message(" Direction: Epigenetic Age Acceleration → Feruli 4 sulfate")
message(strrep("=", 60))

# -------------------- Load Metabolite GWAS as Outcome --------------------
message("\n Loading metabolite GWAS data as outcome...")

# Load full metabolite GWAS data
metabolite_gwas <- fread("/vol/projects/yzhang/500FG_aging/input/MR/ferulic_4_sulfate/phenocode-X100005389_Plasma_NHW.tsv.gz")

message("  ✓ Loaded ", nrow(metabolite_gwas), " SNPs from metabolite GWAS")

# Format as outcome data
metabolite_out_dat <- metabolite_gwas %>%
  as.data.frame() %>%
  format_data(
    type = "outcome",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    pval_col = "pval",
    chr_col = "chrom",
    pos_col = "pos"
  ) %>%
  mutate(outcome = "Feruli 4 sulfate")

if (!"eaf.outcome" %in% names(metabolite_out_dat)) metabolite_out_dat$eaf.outcome <- NA_real_

message("  ✓ Metabolite outcome data formatted: ", nrow(metabolite_out_dat), " SNPs")

# -------------------- Define Reverse MR Parameters --------------------
N_METABOLITE <- 2392  # Metabolite sample size (outcome in reverse MR)
P_THRESHOLD_EXP <- 5e-8  # P-value threshold for exposure instruments

# Initialize results storage for reverse MR
reverse_results <- list()
reverse_mr_res <- list()

# Output directory for reverse MR
reverse_output_dir <- file.path(output_base_dir, "reverse_MR")
dir.create(reverse_output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOOP THROUGH ALL EPIGENETIC CLOCKS AS EXPOSURES (REVERSE MR)
# ============================================================

for (i in seq_along(outcomes_list)) {
  
  clock_info <- outcomes_list[[i]]
  clock_name <- clock_info$name
  clock_file <- clock_info$file
  N_CLOCK <- clock_info$sample_size
  
  message("\n", strrep("=", 60))
  message(" REVERSE MR ", i, "/", length(outcomes_list), ": ", clock_name, " → Feruli 4 sulfate")
  message(strrep("=", 60))
  
  tryCatch({
    
    # -------------------- Load Clock GWAS Data --------------------
    message("\n Loading clock GWAS data: ", clock_file)
    clock_data <- readRDS(clock_file)
    message("  ✓ Loaded ", nrow(clock_data), " rows")
    
    # -------------------- Filter Significant SNPs --------------------
    message(" Filtering significant SNPs (p < ", P_THRESHOLD_EXP, ")...")
    
    clock_sig <- clock_data %>%
      filter(P < P_THRESHOLD_EXP)
    
    n_sig <- nrow(clock_sig)
    message("  ✓ Found ", n_sig, " genome-wide significant SNPs")
    
    if (n_sig == 0) {
      message("  ⚠️ No significant SNPs found for ", clock_name, ". Trying relaxed threshold (p < 5e-6)...")
      
      clock_sig <- clock_data %>%
        filter(P < 5e-6)
      
      n_sig <- nrow(clock_sig)
      message("  ✓ Found ", n_sig, " suggestive SNPs (p < 5e-6)")
      
      if (n_sig == 0) {
        message("   Still no SNPs. Skipping ", clock_name)
        next
      }
    }
    
    # -------------------- Format as Exposure Data --------------------
    message(" Formatting as exposure data...")
    
    clock_exp_dat <- clock_sig %>%
      rename(
        effect_allele = A1,
        other_allele = A2,
        beta = Effect,
        se = SE,
        pval = P
      ) %>%
      format_data(
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        chr_col = "chr",
        pos_col = "pos"
      ) %>%
      mutate(exposure = clock_name)
    
    if (!"eaf.exposure" %in% names(clock_exp_dat)) clock_exp_dat$eaf.exposure <- NA_real_
    
    message("  ✓ Exposure SNPs before clumping: ", nrow(clock_exp_dat))
    
    # -------------------- LD Clumping --------------------
    message(" Performing LD clumping...")
    
    clock_exp_clumped <- tryCatch(
      clump_data(clock_exp_dat, clump_r2 = 0.001, clump_kb = 10000),
      error = function(e) {
        message("  ⚠️ Online clumping failed: ", e$message)
        message("  → Using local LD pruning list...")
        
        # Use local LD-pruned SNPs if online clumping fails
        snps_ind <- read.table("/vol/projects/mballan/send/EUR_pruned.prune.in", 
                               header = FALSE, stringsAsFactors = FALSE)
        clock_exp_dat %>% filter(SNP %in% snps_ind$V1)
      }
    )
    
    n_clumped <- nrow(clock_exp_clumped)
    message("  ✓ Independent SNPs after clumping: ", n_clumped)
    
    if (n_clumped == 0) {
      message("   No independent SNPs after clumping. Skipping ", clock_name)
      next
    }
    
    # -------------------- Harmonise Data --------------------
    message(" Harmonising with metabolite data...")
    
    dat_reverse <- harmonise_data(clock_exp_clumped, metabolite_out_dat, action = 2)
    
    # Handle list columns
    if ("remove" %in% names(dat_reverse) && is.list(dat_reverse$remove)) {
      dat_reverse$remove <- unlist(dat_reverse$remove)
    }
    if ("ambiguous" %in% names(dat_reverse) && is.list(dat_reverse$ambiguous)) {
      dat_reverse$ambiguous <- unlist(dat_reverse$ambiguous)
    }
    
    dat_reverse <- dat_reverse %>%
      filter(remove == FALSE, is.na(ambiguous) | ambiguous == FALSE)
    
    if (nrow(dat_reverse) == 0) {
      message("   No overlapping SNPs after harmonisation. Skipping ", clock_name)
      next
    }
    
    message("  ✓ Harmonised SNPs: ", nrow(dat_reverse))
    
    dat_reverse$samplesize.exposure <- N_CLOCK
    dat_reverse$samplesize.outcome <- N_METABOLITE
    
    # -------------------- Calculate F-statistics --------------------
    message(" Calculating F-statistics...")
    
    dat_reverse <- dat_reverse %>%
      mutate(F_stat = (beta.exposure^2) / (se.exposure^2))
    
    F_mean_rev <- mean(dat_reverse$F_stat, na.rm = TRUE)
    message("  ✓ Mean F-statistic: ", round(F_mean_rev, 2))
    
    # Filter by F-statistic
    dat_reverse_filtered <- dat_reverse %>% filter(F_stat >= F_MIN)
    
    if (nrow(dat_reverse_filtered) == 0) {
      message("   No SNPs with F >= ", F_MIN, ". Skipping ", clock_name)
      next
    }
    
    message("  ✓ SNPs after F-stat filter: ", nrow(dat_reverse_filtered))
    
    # -------------------- Run MR Analyses --------------------
    message(" Running MR analyses...")
    
    ns_rev <- length(unique(dat_reverse_filtered$SNP))
    message("  ✓ Number of SNPs used: ", ns_rev)
    
    # Define methods
    mlist_rev <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median")
    if (has_raps) mlist_rev <- c(mlist_rev, "mr_raps")
    
    # Run MR
    mr_res_rev <- tryCatch(
      mr(dat_reverse_filtered, method_list = mlist_rev),
      error = function(e) {
        message("  ⚠️ MR failed: ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(mr_res_rev)) {
      message("   MR analysis failed. Skipping ", clock_name)
      next
    }
    
    # Extract results
    ivw_rev <- pickrow(mr_res_rev, "Inverse variance weighted")
    wald_rev <- pickrow(mr_res_rev, "^Wald ratio$")
    
    if (ns_rev == 1 && !is.null(wald_rev)) {
      method_rev <- "Wald ratio"
      beta_rev <- safe1(wald_rev, "b")
      se_rev <- safe1(wald_rev, "se")
      p_rev <- safe1(wald_rev, "pval")
    } else {
      method_rev <- "Inverse variance weighted (random effects)"
      beta_rev <- safe1(ivw_rev, "b")
      se_rev <- safe1(ivw_rev, "se")
      p_rev <- safe1(ivw_rev, "pval")
    }
    
    message("  ✓ Reverse MR: Beta = ", round(beta_rev, 4),
            ", SE = ", round(se_rev, 4),
            ", P = ", format(p_rev, scientific = TRUE, digits = 3))
    
    # -------------------- Sensitivity Tests --------------------
    message(" Running sensitivity tests...")
    
    # Heterogeneity
    het_rev <- tryCatch(mr_heterogeneity(dat_reverse_filtered), error = function(e) NULL)
    Q_p_rev <- NA_real_
    if (!is.null(het_rev) && nrow(het_rev)) {
      hi_rev <- het_rev[het_rev$method == "Inverse variance weighted", , drop = FALSE]
      if (nrow(hi_rev)) Q_p_rev <- as.numeric(hi_rev$Q_pval[1])
    }
    
    # Pleiotropy
    pleio_rev <- tryCatch(mr_pleiotropy_test(dat_reverse_filtered), error = function(e) NULL)
    egger_p_rev <- if (!is.null(pleio_rev) && nrow(pleio_rev)) as.numeric(pleio_rev$pval[1]) else NA_real_
    
    # MR-PRESSO
    presso_p_rev <- NA_real_
    if (has_mrpresso && ns_rev >= 3) {
      pr_rev <- tryCatch(
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = as.data.frame(dat_reverse_filtered),
          NbDistribution = 1000,
          SignifThreshold = 0.05
        ),
        error = function(e) NULL
      )
      if (!is.null(pr_rev) && !is.null(pr_rev$`MR-PRESSO results`)) {
        presso_p_rev <- tryCatch(pr_rev$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
      }
    }
    
    # Extract other methods
    raps_rev <- pickrow(mr_res_rev, "RAPS")
    egger_rev <- pickrow(mr_res_rev, "MR Egger")
    wm_rev <- pickrow(mr_res_rev, "Weighted median")
    
    # -------------------- Compile Results --------------------
    results_rev <- tibble::tibble(
      direction = "Reverse",
      exposure = clock_name,
      outcome = "Feruli 4 sulfate",
      n_SNPs_used = ns_rev,
      F_mean = F_mean_rev,
      Q_pval_heterogeneity = Q_p_rev,
      Egger_intercept_pval = egger_p_rev,
      MR_PRESSO_pval = presso_p_rev,
      MR_method = method_rev,
      MR_beta = beta_rev,
      MR_se = se_rev,
      MR_pval = p_rev,
      MR_Egger_beta = safe1(egger_rev, "b"),
      MR_Egger_pval = safe1(egger_rev, "pval"),
      MR_WM_beta = safe1(wm_rev, "b"),
      MR_WM_pval = safe1(wm_rev, "pval"),
      MR_RAPS_beta = safe1(raps_rev, "b"),
      MR_RAPS_pval = safe1(raps_rev, "pval")
    )
    
    # Store results
    reverse_results[[clock_name]] <- results_rev
    reverse_mr_res[[clock_name]] <- list(
      mr_res = mr_res_rev,
      het = het_rev,
      pleio = pleio_rev,
      dat_filtered = dat_reverse_filtered
    )
    
    # -------------------- Save Individual Results --------------------
    rev_out_dir <- file.path(reverse_output_dir, paste0("reverse_MR_", clock_name))
    dir.create(rev_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    write.csv(results_rev, file.path(rev_out_dir, "reverse_MR_results.csv"), row.names = FALSE)
    write.csv(mr_res_rev, file.path(rev_out_dir, "reverse_MR_all_methods.csv"), row.names = FALSE)
    write.csv(dat_reverse_filtered, file.path(rev_out_dir, "reverse_harmonised_data.csv"), row.names = FALSE)
    
    # -------------------- Generate Plots --------------------
    message(" Generating plots...")
    
    tryCatch({
      p_scatter <- mr_scatter_plot(mr_res_rev, dat_reverse_filtered)
      ggsave(file.path(rev_out_dir, "reverse_scatter_plot.png"), p_scatter[[1]],
             width = 10, height = 8, dpi = 300)
    }, error = function(e) message("  ⚠️ Scatter plot failed"))
    
    tryCatch({
      single_snp <- mr_singlesnp(dat_reverse_filtered)
      p_forest <- mr_forest_plot(single_snp)
      ggsave(file.path(rev_out_dir, "reverse_forest_plot.png"), p_forest[[1]],
             width = 10, height = 10, dpi = 300)
    }, error = function(e) message("  ⚠️ Forest plot failed"))
    
    message("  ✓ Completed reverse MR for: ", clock_name)
    
  }, error = function(e) {
    message("   ERROR processing ", clock_name, ": ", e$message)
  })
}

message("\n", strrep("=", 60))
message(" REVERSE MR LOOP COMPLETED")
message(strrep("=", 60))

# ============================================================
# BIDIRECTIONAL MR VISUALIZATION
# ============================================================

message(" Creating bidirectional comparison plot...")

tryCatch({
  
  if (exists("bidirectional_results") && nrow(bidirectional_results) > 0) {
    
    # Prepare data for plotting
    plot_data <- bidirectional_results %>%
      mutate(
        label = paste0(direction, ": ", exposure, " → ", outcome),
        direction = factor(direction, levels = c("Forward", "Reverse")),
        pval_label = ifelse(pval < 0.001, 
                            format(pval, scientific = TRUE, digits = 2),
                            round(pval, 3)),
        sig_color = case_when(
          pval < 0.05/8 ~ "Bonferroni significant",
          pval < 0.05 ~ "Nominally significant", 
          TRUE ~ "Not significant"
        )
      )
    
    # Create bidirectional forest plot
    p_bidirectional <- ggplot(plot_data, aes(x = beta, y = reorder(label, -as.numeric(direction)))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.25, size = 0.8) +
      geom_point(aes(color = sig_color, shape = direction), size = 4) +
      geom_text(aes(label = paste0("p=", pval_label)), 
                hjust = -0.3, vjust = -0.5, size = 3) +
      scale_color_manual(
        values = c("Bonferroni significant" = "#E63946",
                   "Nominally significant" = "#F4A261",
                   "Not significant" = "gray50"),
        name = "Significance"
      ) +
      scale_shape_manual(
        values = c("Forward" = 16, "Reverse" = 17),
        name = "Direction"
      ) +
      labs(
        title = "Bidirectional Mendelian Randomization Results",
        subtitle = "Forward: Ferulic acid → Aging | Reverse: Aging → Ferulic acid",
        x = "MR Effect Estimate (Beta) with 95% CI",
        y = ""
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10)
      ) +
      facet_wrap(~ direction, scales = "free_y", ncol = 1)
    
    # Save plot
    ggsave(
      file.path(output_base_dir, "bidirectional_MR_forest_plot.pdf"),
      p_bidirectional,
      width = 12, height = 10, bg = "white"
    )
    
    ggsave(
      file.path(output_base_dir, "bidirectional_MR_forest_plot.png"),
      p_bidirectional,
      width = 12, height = 10, dpi = 300, bg = "white"
    )
    
    message(" ✓ Bidirectional forest plot saved")
    
    # Display plot
    print(p_bidirectional)
    
  } else {
    message(" ⚠️ No bidirectional results available for plotting")
  }
  
}, error = function(e) {
  message("  Plot generation failed: ", e$message)
})

message("\n ALL BIDIRECTIONAL MR ANALYSES COMPLETE!")

# ============================================================
# BIDIRECTIONAL COMPARISON
# ============================================================

# 整理正向结果
forward_summary <- combined_results %>%
  transmute(direction = "Forward", exposure = "Feruli 4 sulfate", outcome = outcome,
            n_SNPs = n_SNPs_used, beta = MR_beta_estimate, se = MR_standard_error, pval = MR_pval_main)

# 合并双向结果
if (exists("combined_reverse") && nrow(combined_reverse) > 0) {
  bidirectional <- bind_rows(forward_summary, combined_reverse)
  
  cat("\n========== BIDIRECTIONAL MR COMPARISON ==========\n")
  print(bidirectional %>% mutate(pval = format(pval, scientific = TRUE, digits = 2)))
  
  # 保存
  write.csv(bidirectional, file.path(output_base_dir, "bidirectional_comparison.csv"), row.names = FALSE)
  
  # 解读
  fwd_sig <- sum(forward_summary$pval < 0.05)
  rev_sig <- sum(combined_reverse$pval < 0.05)
  
  cat("\n========== INTERPRETATION ==========\n")
  if (fwd_sig > 0 && rev_sig == 0) {
    cat("✓ Forward significant, Reverse not → Supports Ferulic acid → Aging causality\n")
  } else if (fwd_sig == 0 && rev_sig > 0) {
    cat("⚠ Reverse significant only → Possible reverse causation\n")
  } else if (fwd_sig > 0 && rev_sig > 0) {
    cat("⚠ Both directions significant → Possible bidirectional causation\n")
  } else {
    cat("○ Neither direction significant → No causal evidence\n")
  }
} else {
  cat("\nNo reverse MR results - epigenetic clocks likely lack GWAS-significant SNPs\n")
  cat("This itself suggests reverse causation is unlikely (no strong genetic instruments for aging)\n")
}

save.image("/vol/projects/yzhang/500FG_aging/output/16_MR/ferulic_4_sulfate/MR_ferulic_4_sulfate_reverse.Rdata")
